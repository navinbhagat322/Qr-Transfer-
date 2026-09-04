#!/usr/bin/env python3
"""Pure stdlib QR decoder for the offline_qr_gen.py "PART i/N:<base64>" transfer format.

No third-party packages are used anywhere in this file (only: sys, os, glob, zlib,
struct, argparse, base64, math). Verify with `grep -n '^import\|^from' this file.

Supported image formats
------------------------
  * PPM (P6, binary RGB) -- exactly what offline_qr_gen.py's render_ppm() writes.
  * PNG (8-bit/channel, color types 0 gray / 2 RGB / 3 palette / 6 RGBA,
    non-interlaced). Decoded with a hand-rolled chunk/filter parser; the
    compressed IDAT stream itself is inflated with the stdlib `zlib` module
    (zlib/DEFLATE decompression is not "a QR/image codec", it's generic
    stdlib compression, so this stays dependency-free).
  * PBM (P4, binary bitmap) and PGM (P5, binary grayscale) are supported
    trivially as a bonus since they share PPM's binary-image plumbing.

Not supported / known limitations
----------------------------------
  * PNG interlacing (Adam7) is NOT implemented -- interlaced PNGs will raise.
  * 16-bit-per-channel PNG, colortype 4 (gray+alpha) are NOT implemented.
  * Only axis-aligned, non-rotated, non-skewed, reasonably sharp images are
    supported. There is NO perspective correction and no deskew step. A
    photo taken at an angle, or a heavily blurred/rotated capture, will very
    likely fail to decode. This matches the "clean, well-aligned rendered
    image" requirement; robust real-world camera capture is out of scope.
  * Reed-Solomon error correction is FULLY implemented (syndromes +
    Berlekamp-Massey + Chien search + Forney algorithm), adapted from the
    same GF(256) arithmetic offline_qr_gen.py uses for encoding. It can
    correct up to floor(ecc_codewords/2) byte errors per block, exactly as
    the QR spec intends -- this is not merely a checksum/verify path.
  * Only QR "byte mode" symbols are parsed (that's all offline_qr_gen.py
    ever produces). Other encoding modes (numeric/alphanumeric/kanji) are
    not decoded.
  * Version is derived purely from the measured symbol size (distance
    between finder patterns), not from reading/BCH-decoding the version
    info bits -- reliable for the axis-aligned case this tool targets.
    Format info (EC level + mask) IS decoded from the image, using
    brute-force matching against all 32 possible 15-bit codewords (which
    is an exact substitute for BCH(15,5) decoding since the codeword space
    is tiny).

CLI usage
---------
    python3 laptop_qr_decoder.py --input-dir ppm_out --output reconstructed.bin

    python3 laptop_qr_decoder.py -i /path/to/frames -o out/myfile.zip [--verbose]

The script reads every .ppm/.pgm/.pbm/.png file in --input-dir, decodes each
as a QR symbol, extracts "PART i/N:<b64chunk>" payloads, sorts by i, verifies
all N parts of a single "generation" are present, concatenates the base64
text in order, base64-decodes it, and writes the resulting bytes to --output.
"""
import sys
import os
import glob
import zlib
import struct
import base64
import argparse

# =========================================================================
# GF(256) arithmetic -- identical field/generator to offline_qr_gen.py
# =========================================================================

_EXP = [0] * 512
_LOG = [0] * 256


def _init_gf():
    x = 1
    for i in range(255):
        _EXP[i] = x
        _LOG[x] = i
        x <<= 1
        if x & 0x100:
            x ^= 0x11D
    for i in range(255, 512):
        _EXP[i] = _EXP[i - 255]


_init_gf()


def gf_mul(a, b):
    if a == 0 or b == 0:
        return 0
    return _EXP[_LOG[a] + _LOG[b]]


def gf_pow(a, n):
    if n == 0:
        return 1
    return _EXP[(_LOG[a] * n) % 255]


def gf_inv(a):
    return _EXP[255 - _LOG[a]]


def gf_poly_scale(p, x):
    return [gf_mul(c, x) for c in p]


def gf_poly_add(p, q):
    r = [0] * max(len(p), len(q))
    for i, c in enumerate(p):
        r[i + len(r) - len(p)] ^= c
    for i, c in enumerate(q):
        r[i + len(r) - len(q)] ^= c
    return r


def gf_poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, pc in enumerate(p):
        if pc == 0:
            continue
        for j, qc in enumerate(q):
            r[i + j] ^= gf_mul(pc, qc)
    return r


def gf_poly_eval(p, x):
    y = p[0]
    for c in p[1:]:
        y = gf_mul(y, x) ^ c
    return y


def rs_generator_poly(n_ecc):
    g = [1]
    for i in range(n_ecc):
        g.append(0)
        for j in range(len(g) - 1, 0, -1):
            g[j] ^= gf_mul(g[j - 1], _EXP[i])
    return g


# ------------------------- RS syndrome decoding ---------------------------
# Full error-correcting decode: syndromes -> Berlekamp-Massey -> Chien search
# -> Forney algorithm. Roots of the generator are alpha^0..alpha^(nsym-1),
# matching offline_qr_gen.py's rs_generator_poly exactly.

class RSDecodeError(Exception):
    pass


def rs_calc_syndromes(msg, nsym):
    # msg is coefficients, highest degree first (msg[0] = most significant byte)
    synd = [0] * nsym
    for i in range(nsym):
        synd[i] = gf_poly_eval(msg, gf_pow(2, i))
    return synd


def rs_find_error_locator(synd):
    err_loc = [1]
    old_loc = [1]
    for i in range(len(synd)):
        old_loc = old_loc + [0]
        delta = synd[i]
        for j in range(1, len(err_loc)):
            delta ^= gf_mul(err_loc[-(j + 1)], synd[i - j])
        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = gf_poly_scale(old_loc, delta)
                old_loc = gf_poly_scale(err_loc, gf_inv(delta))
                err_loc = new_loc
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))
    # strip leading zeros
    while err_loc and err_loc[0] == 0:
        err_loc.pop(0)
    errs = len(err_loc) - 1
    return err_loc, errs


def rs_find_errors(err_loc, n):
    errs = len(err_loc) - 1
    err_pos = []
    for i in range(n):
        if gf_poly_eval(err_loc, gf_pow(2, i)) == 0:
            err_pos.append(n - 1 - i)
    if len(err_pos) != errs:
        raise RSDecodeError("could not locate all errors (uncorrectable block)")
    return err_pos


def rs_correct_errata(msg, synd, err_pos):
    n = len(msg)
    coef_pos = [n - 1 - p for p in err_pos]
    err_loc = [1]
    for p in coef_pos:
        err_loc = gf_poly_mul(err_loc, [gf_pow(2, p), 1])
    # error evaluator polynomial: synd * err_loc mod x^nsym, then reversed
    synd_rev = synd[::-1]
    err_eval = gf_poly_mul(synd_rev, err_loc)
    err_eval = err_eval[len(err_eval) - len(synd):]

    err_loc_prime_tmp = []
    for i in range(len(coef_pos)):
        x = gf_pow(2, coef_pos[i])
        e = 1
        for j in range(len(coef_pos)):
            if j != i:
                e = gf_mul(e, 1 ^ gf_mul(x, gf_pow(2, coef_pos[j])))
        err_loc_prime_tmp.append(e)

    x_list = [gf_pow(2, p) for p in coef_pos]
    e = [0] * n
    for i, xv in enumerate(x_list):
        xv_inv = gf_inv(xv)
        y = gf_poly_eval(err_eval[::-1], xv_inv)
        y = gf_mul(gf_pow(xv, 1), y)
        magnitude = gf_mul(y, gf_inv(err_loc_prime_tmp[i]))
        e[coef_pos[i]] = magnitude
    msg = [c ^ e[n - 1 - i] if False else c for i, c in enumerate(msg)]
    out = list(msg)
    for p in err_pos:
        idx = p  # position from start, 0 = first byte
        coef_index = n - 1 - idx
        out[idx] ^= e[coef_index]
    return out


def rs_decode(codewords, nsym):
    """codewords: full block (data+ecc) as list of ints, data first (matches
    offline_qr_gen.py's block layout: data codewords followed by ecc codewords).
    Returns corrected data codewords (without ecc) or raises RSDecodeError."""
    msg = list(codewords)
    synd = rs_calc_syndromes(msg, nsym)
    if max(synd) == 0:
        return msg[:-nsym]  # no errors
    err_loc, errs = rs_find_error_locator(synd)
    if errs * 2 > nsym:
        raise RSDecodeError("too many errors to correct in this block")
    err_pos = rs_find_errors(err_loc, len(msg))
    corrected = rs_correct_errata(msg, synd, err_pos)
    # verify
    if max(rs_calc_syndromes(corrected, nsym)) != 0:
        raise RSDecodeError("RS correction failed verification")
    return corrected[:-nsym]


# =========================================================================
# QR structural tables (copied from offline_qr_gen.py so this file is
# self-contained / standalone)
# =========================================================================

QR_L = {
    1: (26, 7, 1, 19, 0, 0), 2: (44, 10, 1, 34, 0, 0), 3: (70, 15, 1, 55, 0, 0),
    4: (100, 20, 1, 80, 0, 0), 5: (134, 26, 1, 108, 0, 0), 6: (172, 18, 2, 68, 0, 0),
    7: (196, 20, 2, 78, 0, 0), 8: (242, 24, 2, 97, 0, 0), 9: (292, 30, 2, 116, 0, 0),
    10: (346, 18, 2, 68, 2, 69), 11: (404, 20, 4, 81, 0, 0), 12: (466, 24, 2, 92, 2, 93),
    13: (532, 26, 4, 107, 0, 0), 14: (581, 30, 3, 115, 1, 116), 15: (655, 22, 5, 87, 1, 88),
    16: (733, 24, 5, 98, 1, 99), 17: (815, 28, 1, 107, 5, 108), 18: (901, 30, 5, 120, 1, 121),
    19: (991, 28, 3, 113, 4, 114), 20: (1085, 28, 3, 107, 5, 108), 21: (1156, 28, 4, 116, 4, 117),
    22: (1258, 28, 2, 111, 7, 112), 23: (1364, 30, 4, 121, 5, 122), 24: (1474, 30, 6, 117, 4, 118),
    25: (1588, 26, 8, 106, 4, 107), 26: (1706, 28, 10, 114, 2, 115), 27: (1828, 30, 8, 122, 4, 123),
    28: (1921, 30, 3, 117, 10, 118), 29: (2051, 30, 7, 116, 7, 117), 30: (2185, 30, 5, 115, 10, 116),
    31: (2323, 30, 13, 115, 3, 116), 32: (2465, 30, 17, 115, 0, 0), 33: (2611, 30, 17, 115, 1, 116),
    34: (2761, 30, 13, 115, 6, 116), 35: (2876, 30, 12, 121, 7, 122), 36: (2996, 30, 6, 121, 14, 122),
    37: (3122, 30, 17, 122, 4, 123), 38: (3251, 30, 4, 122, 18, 123), 39: (3388, 30, 20, 117, 4, 118),
    40: (3417, 30, 19, 118, 6, 119),
}

ALIGNMENT_COORDS = {
    1: [], 2: [6, 18], 3: [6, 22], 4: [6, 26], 5: [6, 30], 6: [6, 34],
    7: [6, 22, 38], 8: [6, 24, 42], 9: [6, 26, 46], 10: [6, 28, 50],
    11: [6, 30, 54], 12: [6, 32, 58], 13: [6, 34, 62], 14: [6, 26, 46, 66],
    15: [6, 26, 48, 70], 16: [6, 26, 50, 74], 17: [6, 30, 54, 78],
    18: [6, 30, 56, 82], 19: [6, 30, 58, 86], 20: [6, 34, 62, 90],
    21: [6, 28, 50, 72, 94], 22: [6, 26, 50, 74, 98], 23: [6, 30, 54, 78, 102],
    24: [6, 28, 54, 80, 106], 25: [6, 32, 58, 84, 110], 26: [6, 30, 58, 86, 114],
    27: [6, 34, 62, 90, 118], 28: [6, 26, 50, 74, 98, 122], 29: [6, 30, 54, 78, 102, 126],
    30: [6, 26, 52, 78, 104, 130], 31: [6, 30, 56, 82, 108, 134], 32: [6, 34, 60, 86, 112, 138],
    33: [6, 30, 58, 86, 114, 142], 34: [6, 34, 62, 90, 118, 146], 35: [6, 30, 54, 78, 102, 126, 150],
    36: [6, 24, 50, 76, 102, 128, 154], 37: [6, 28, 54, 80, 106, 132, 158],
    38: [6, 32, 58, 84, 110, 136, 162], 39: [6, 26, 54, 82, 110, 138, 166],
    40: [6, 30, 58, 86, 114, 142, 170],
}

REMAINDER_BITS = [0, 0, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3,
                  4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0]

EC_LEVEL_INDICATOR = {'L': 0b01, 'M': 0b00, 'Q': 0b11, 'H': 0b10}
EC_LEVEL_FROM_INDICATOR = {v: k for k, v in EC_LEVEL_INDICATOR.items()}

FORMAT_MASK = 0b101010000010010


def char_count_bits(version):
    return 8 if version <= 9 else 16


def bch_format_info(fmt5):
    g = 0b10100110111
    data = fmt5 << 10
    val = data
    while val.bit_length() > 10:
        val ^= g << (val.bit_length() - 11)
    return (data | val) ^ FORMAT_MASK


def mask_condition(mask_id, r, c):
    if mask_id == 0:
        return (r + c) % 2 == 0
    if mask_id == 1:
        return r % 2 == 0
    if mask_id == 2:
        return c % 3 == 0
    if mask_id == 3:
        return (r + c) % 3 == 0
    if mask_id == 4:
        return (r // 2 + c // 3) % 2 == 0
    if mask_id == 5:
        return (r * c) % 2 + (r * c) % 3 == 0
    if mask_id == 6:
        return ((r * c) % 2 + (r * c) % 3) % 2 == 0
    if mask_id == 7:
        return ((r + c) % 2 + (r * c) % 3) % 2 == 0
    raise ValueError("bad mask id")


class QRSkeleton:
    """Rebuilds the function-module layout (finder/timing/alignment/format
    reservations) for a given version, mirroring offline_qr_gen.py's
    QRMatrix builder minus the actual data placement. Used to know which
    modules are function modules (skip when reading data) and to walk the
    same zig-zag data order the encoder used, in reverse (i.e. for reading).
    """

    def __init__(self, version):
        self.version = version
        self.size = version * 4 + 17
        self.is_function = [[False] * self.size for _ in range(self.size)]

    def mark(self, r, c):
        if 0 <= r < self.size and 0 <= c < self.size:
            self.is_function[r][c] = True

    def place_finder(self, r, c):
        for dr in range(-1, 8):
            for dc in range(-1, 8):
                rr, cc = r + dr, c + dc
                if 0 <= rr < self.size and 0 <= cc < self.size:
                    self.mark(rr, cc)

    def build(self):
        self.place_finder(0, 0)
        self.place_finder(0, self.size - 7)
        self.place_finder(self.size - 7, 0)
        coords = ALIGNMENT_COORDS[self.version]
        for r in coords:
            for c in coords:
                if (r <= 8 and c <= 8) or (r <= 8 and c >= self.size - 9) or (r >= self.size - 9 and c <= 8):
                    continue
                for dr in range(-2, 3):
                    for dc in range(-2, 3):
                        self.mark(r + dr, c + dc)
        for i in range(8, self.size - 8):
            self.mark(6, i)
            self.mark(i, 6)
        self.mark(self.size - 8, 8)  # dark module
        for i in range(9):
            self.mark(8, i)
            self.mark(i, 8)
        for i in range(self.size - 8, self.size):
            self.mark(8, i)
        for i in range(self.size - 7, self.size):
            self.mark(i, 8)
        if self.version >= 7:
            for r in range(6):
                for c in range(self.size - 11, self.size - 8):
                    self.mark(r, c)
            for c in range(6):
                for r in range(self.size - 11, self.size - 8):
                    self.mark(r, c)
        return self

    def data_order(self):
        """Return list of (row,col) in the exact order the encoder's
        place_data() consumes bits -- i.e. the order to read data bits."""
        order = []
        col = self.size - 1
        upward = True
        while col > 0:
            if col == 6:
                col -= 1
            for i in range(self.size):
                row = (self.size - 1 - i) if upward else i
                for c in (col, col - 1):
                    if not self.is_function[row][c]:
                        order.append((row, c))
            upward = not upward
            col -= 2
        return order


# =========================================================================
# Image loading: PPM/PGM/PBM and PNG, stdlib only
# =========================================================================

class Image:
    """Minimal RGB image wrapper: width, height, get(x,y) -> (r,g,b)."""

    def __init__(self, width, height, pixels):
        self.width = width
        self.height = height
        self.pixels = pixels  # flat bytearray, 3 bytes per pixel (RGB)

    def get_gray(self, x, y):
        i = (y * self.width + x) * 3
        p = self.pixels
        return (p[i] * 299 + p[i + 1] * 587 + p[i + 2] * 114) // 1000


def _read_pnm_header(f):
    """Read a PNM header (magic, width, height, maxval[for P5/P6]),
    tolerating whitespace and '#' comments per the NetPBM spec."""
    def token():
        buf = bytearray()
        while True:
            ch = f.read(1)
            if not ch:
                raise ValueError("unexpected EOF in PNM header")
            if ch in b'#':
                while ch and ch != b'\n':
                    ch = f.read(1)
                continue
            if ch.isspace():
                if buf:
                    return bytes(buf)
                continue
            buf += ch

    magic = token().decode('ascii')
    width = int(token())
    height = int(token())
    maxval = None
    if magic in ('P5', 'P6'):
        maxval = int(token())
    return magic, width, height, maxval


def load_pnm(path):
    with open(path, 'rb') as f:
        magic, width, height, maxval = _read_pnm_header(f)
        raw = f.read()
    pixels = bytearray(width * height * 3)
    if magic == 'P6':
        # binary RGB, assume maxval 255 (offline_qr_gen.py always uses 255)
        for i in range(width * height):
            pixels[i * 3:i * 3 + 3] = raw[i * 3:i * 3 + 3]
    elif magic == 'P5':
        for i in range(width * height):
            v = raw[i]
            pixels[i * 3:i * 3 + 3] = bytes((v, v, v))
    elif magic == 'P4':
        # packed 1-bit-per-pixel, MSB first, 1 = black
        row_bytes = (width + 7) // 8
        for y in range(height):
            for x in range(width):
                byte = raw[y * row_bytes + x // 8]
                bit = (byte >> (7 - (x % 8))) & 1
                v = 0 if bit else 255
                idx = (y * width + x) * 3
                pixels[idx:idx + 3] = bytes((v, v, v))
    else:
        raise ValueError(f"unsupported PNM magic {magic!r}")
    return Image(width, height, pixels)


_PAETH = None


def _paeth_predictor(a, b, c):
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def load_png(path):
    with open(path, 'rb') as f:
        sig = f.read(8)
        if sig != b'\x89PNG\r\n\x1a\n':
            raise ValueError("not a PNG file")
        width = height = bitdepth = colortype = None
        interlace = 0
        idat = bytearray()
        palette = None
        while True:
            hdr = f.read(8)
            if len(hdr) < 8:
                break
            length, ctype = struct.unpack('>I4s', hdr)
            data = f.read(length)
            f.read(4)  # crc, unchecked
            ctype = ctype.decode('ascii')
            if ctype == 'IHDR':
                (width, height, bitdepth, colortype, _comp, _filt, interlace) = \
                    struct.unpack('>IIBBBBB', data)
            elif ctype == 'PLTE':
                palette = [tuple(data[i:i + 3]) for i in range(0, len(data), 3)]
            elif ctype == 'IDAT':
                idat += data
            elif ctype == 'IEND':
                break
    if bitdepth != 8:
        raise ValueError(f"unsupported PNG bit depth {bitdepth} (only 8-bit supported)")
    if interlace != 0:
        raise ValueError("interlaced PNG not supported")
    if colortype not in (0, 2, 3, 6):
        raise ValueError(f"unsupported PNG color type {colortype}")

    channels = {0: 1, 2: 3, 3: 1, 6: 4}[colortype]
    raw = zlib.decompress(bytes(idat))
    stride = width * channels
    out = bytearray(width * height * 3)
    prev_row = bytearray(stride)
    pos = 0
    for y in range(height):
        filt = raw[pos]
        pos += 1
        row = bytearray(raw[pos:pos + stride])
        pos += stride
        for x in range(stride):
            a = row[x - channels] if x >= channels else 0
            b = prev_row[x]
            c = prev_row[x - channels] if x >= channels else 0
            if filt == 0:
                pass
            elif filt == 1:
                row[x] = (row[x] + a) & 0xFF
            elif filt == 2:
                row[x] = (row[x] + b) & 0xFF
            elif filt == 3:
                row[x] = (row[x] + (a + b) // 2) & 0xFF
            elif filt == 4:
                row[x] = (row[x] + _paeth_predictor(a, b, c)) & 0xFF
            else:
                raise ValueError(f"unsupported PNG filter type {filt}")
        # convert this row to RGB into out
        for x in range(width):
            base = x * channels
            if colortype == 0:
                v = row[base]
                rgb = (v, v, v)
            elif colortype == 2:
                rgb = (row[base], row[base + 1], row[base + 2])
            elif colortype == 3:
                rgb = palette[row[base]]
            elif colortype == 6:
                rgb = (row[base], row[base + 1], row[base + 2])
            oi = (y * width + x) * 3
            out[oi:oi + 3] = bytes(rgb)
        prev_row = row
    return Image(width, height, out)


def load_image(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == '.png':
        return load_png(path)
    if ext in ('.ppm', '.pgm', '.pbm'):
        return load_pnm(path)
    # fall back: sniff signature
    with open(path, 'rb') as f:
        head = f.read(8)
    if head.startswith(b'\x89PNG'):
        return load_png(path)
    if head[:2] in (b'P6', b'P5', b'P4'):
        return load_pnm(path)
    raise ValueError(f"unrecognized image format: {path}")


# =========================================================================
# Binarization
# =========================================================================

def binarize(img):
    """Otsu's method global threshold -> returns a function is_dark(x,y)."""
    hist = [0] * 256
    w, h = img.width, img.height
    # sample a subset of pixels for speed on large images
    step_x = max(1, w // 400)
    step_y = max(1, h // 400)
    total = 0
    for y in range(0, h, step_y):
        for x in range(0, w, step_x):
            hist[img.get_gray(x, y)] += 1
            total += 1
    if total == 0:
        total = 1
    sum_all = sum(i * hist[i] for i in range(256))
    sum_b = 0
    w_b = 0
    max_var = -1
    threshold = 128
    for t in range(256):
        w_b += hist[t]
        if w_b == 0:
            continue
        w_f = total - w_b
        if w_f == 0:
            break
        sum_b += t * hist[t]
        m_b = sum_b / w_b
        m_f = (sum_all - sum_b) / w_f
        var_between = w_b * w_f * (m_b - m_f) ** 2
        if var_between > max_var:
            max_var = var_between
            threshold = t

    def is_dark(x, y):
        if x < 0 or y < 0 or x >= w or y >= h:
            return False
        return img.get_gray(x, y) <= threshold

    return is_dark


# =========================================================================
# Finder pattern detection (axis-aligned only -- see module docstring)
# =========================================================================

def _line_runs(length, sample_fn):
    """sample_fn(i) -> bool (dark?). Returns list of (color, run_length)."""
    runs = []
    cur = sample_fn(0)
    count = 1
    for i in range(1, length):
        v = sample_fn(i)
        if v == cur:
            count += 1
        else:
            runs.append((cur, count))
            cur = v
            count = 1
    runs.append((cur, count))
    return runs


def _find_pattern_centers_in_runs(runs, start_offsets):
    """Scan consecutive runs [dark, light, dark, light, dark] with
    ratio ~1:1:3:1:1 and return list of (center_offset_along_line, module_size)."""
    results = []
    for i in range(len(runs) - 4):
        colors = [runs[i + k][0] for k in range(5)]
        if colors != [True, False, True, False, True]:
            continue
        lens = [runs[i + k][1] for k in range(5)]
        total = sum(lens)
        unit = total / 7.0
        if unit < 1:
            continue
        ok = True
        for k in (0, 1, 3, 4):
            if not (0.5 * unit <= lens[k] <= 1.5 * unit):
                ok = False
                break
        if ok and not (2.5 * unit <= lens[2] <= 3.5 * unit):
            ok = False
        if not ok:
            continue
        offset_start = start_offsets[i]
        center = offset_start + lens[0] + lens[1] + lens[2] / 2.0
        results.append((center, unit))
    return results


def find_finder_centers(is_dark, width, height):
    row_candidates = []  # (x, y, unit)
    for y in range(height):
        runs = _line_runs(width, lambda x, y=y: is_dark(x, y))
        offsets = []
        pos = 0
        for _, ln in runs:
            offsets.append(pos)
            pos += ln
        for cx, unit in _find_pattern_centers_in_runs(runs, offsets):
            row_candidates.append((cx, y, unit))

    # verify each row candidate vertically, and refine y
    confirmed = []
    for cx, cy, unit in row_candidates:
        x = int(round(cx))
        if x < 0 or x >= width:
            continue
        runs = _line_runs(height, lambda y, x=x: is_dark(x, y))
        offsets = []
        pos = 0
        for _, ln in runs:
            offsets.append(pos)
            pos += ln
        col_hits = _find_pattern_centers_in_runs(runs, offsets)
        best = None
        for cy2, unit2 in col_hits:
            if abs(cy2 - cy) <= unit * 3:
                if best is None or abs(cy2 - cy) < abs(best[0] - cy):
                    best = (cy2, unit2)
        if best:
            confirmed.append((cx, best[0], (unit + best[1]) / 2.0))

    # cluster nearby confirmed points
    clusters = []
    for x, y, u in confirmed:
        placed = False
        for cl in clusters:
            if abs(cl['x'] / cl['n'] - x) < u * 2 and abs(cl['y'] / cl['n'] - y) < u * 2:
                cl['x'] += x
                cl['y'] += y
                cl['u'] += u
                cl['n'] += 1
                placed = True
                break
        if not placed:
            clusters.append({'x': x, 'y': y, 'u': u, 'n': 1})

    centers = [(cl['x'] / cl['n'], cl['y'] / cl['n'], cl['u'] / cl['n']) for cl in clusters]
    return centers


def locate_finders(centers):
    """Given >=3 finder-center candidates, pick the top-left, top-right,
    bottom-left triple (axis-aligned assumption)."""
    if len(centers) < 3:
        raise ValueError(f"found only {len(centers)} finder pattern(s), need 3")
    best = _best_candidate_triple(centers)
    if best is None:
        raise ValueError("could not identify finder pattern triple geometry")
    return best[1], best[2], best[3]


def _best_candidate_triple(centers):
    """Return (score, tl, tr, bl) for the single best right-angle triple
    among `centers`, or None if none found. Lower score = more confidently
    a real QR corner triple (close to a right angle AND legs of similar
    length, since a QR's TL-TR and TL-BL finder distances are equal)."""
    import itertools
    best = None
    for combo in itertools.combinations(centers, 3):
        for tl, other1, other2 in itertools.permutations(combo):
            if other1[0] >= tl[0] and other2[1] >= tl[1]:
                tr, bl = other1, other2
            elif other2[0] >= tl[0] and other1[1] >= tl[1]:
                tr, bl = other2, other1
            else:
                continue
            dx1, dy1 = tr[0] - tl[0], tr[1] - tl[1]
            dx2, dy2 = bl[0] - tl[0], bl[1] - tl[1]
            dot = dx1 * dx2 + dy1 * dy2
            len1 = (dx1 ** 2 + dy1 ** 2) ** 0.5
            len2 = (dx2 ** 2 + dy2 ** 2) ** 0.5
            if len1 < 1 or len2 < 1:
                continue
            cos_angle = abs(dot) / (len1 * len2)
            length_ratio = abs(len1 - len2) / max(len1, len2)
            score = cos_angle * 5 + length_ratio
            if best is None or score < best[0]:
                best = (score, tl, tr, bl)
    return best


def locate_all_finder_triples(centers, max_score=0.6):
    """Find as many non-overlapping (tl, tr, bl) finder triples as possible
    among `centers`, for a photo that may contain multiple QR codes tiled
    in a grid. Greedily picks the best-scoring triple, removes its three
    points, and repeats. `max_score` rejects implausible groupings (e.g.
    mixing corners from two different QR codes)."""
    remaining = list(centers)
    triples = []
    while len(remaining) >= 3:
        best = _best_candidate_triple(remaining)
        if best is None or best[0] > max_score:
            break
        _, tl, tr, bl = best
        triples.append((tl, tr, bl))
        used = {tl, tr, bl}
        remaining = [c for c in remaining if c not in used]
    return triples


# =========================================================================
# Full QR symbol decode from a loaded Image
# =========================================================================

class QRDecodeError(Exception):
    pass


def decode_qr_image(img, verbose=False):
    """Decode a photo assumed to contain exactly one QR symbol (original,
    backward-compatible entry point)."""
    is_dark = binarize(img)
    centers = find_finder_centers(is_dark, img.width, img.height)
    tl, tr, bl = locate_finders(centers)
    return decode_qr_from_triple(img, is_dark, tl, tr, bl, verbose=verbose)


def decode_qr_image_multi(img, verbose=False):
    """Decode a photo that may contain MULTIPLE QR symbols tiled together
    (e.g. a grid image from offline_qr_gen.py --grid). Detects all finder
    center candidates once, greedily groups them into independent (tl, tr,
    bl) triples, and decodes each one. Returns a list of decoded byte
    strings (one per QR successfully decoded); QR codes that fail to
    decode are skipped with a printed warning when verbose.

    For a photo with only one QR code, this returns the same single result
    as decode_qr_image() (list of length 1) -- so single-QR photos keep
    working exactly as before via either entry point."""
    is_dark = binarize(img)
    centers = find_finder_centers(is_dark, img.width, img.height)
    triples = locate_all_finder_triples(centers)
    if not triples:
        raise QRDecodeError(f"found only {len(centers)} usable finder pattern(s), could not form any QR triple")
    results = []
    for tl, tr, bl in triples:
        try:
            data = decode_qr_from_triple(img, is_dark, tl, tr, bl, verbose=verbose)
            results.append(data)
        except QRDecodeError as e:
            if verbose:
                print(f"    [skip one grid cell] {e}")
            continue
    if not results:
        raise QRDecodeError("found finder triples but none decoded successfully")
    return results


def decode_qr_from_triple(img, is_dark, tl, tr, bl, verbose=False):
    """Decode a single QR symbol given its (tl, tr, bl) finder centers,
    reusing a precomputed `is_dark` sampler for the whole image (so this
    works whether the image contains one QR or many)."""
    module_size_x = abs(tr[0] - tl[0]) / max(1, round(abs(tr[0] - tl[0]) / ((tl[2] + tr[2]) / 2)))
    # estimate size (modules) from finder center distance: centers are 3.5
    # modules in from each edge, and (size-7) modules apart center-to-center.
    unit = (tl[2] + tr[2] + bl[2]) / 3.0
    dist_x = abs(tr[0] - tl[0])
    dist_y = abs(bl[1] - tl[1])
    size_est_x = round(dist_x / unit) + 7
    size_est_y = round(dist_y / unit) + 7
    size = size_est_x if abs(size_est_x - size_est_y) <= 2 else max(size_est_x, size_est_y)
    if (size - 17) % 4 != 0:
        # snap to nearest valid QR size
        size = 17 + round((size - 17) / 4.0) * 4
    version = (size - 17) // 4
    if not (1 <= version <= 40):
        raise QRDecodeError(f"implausible version derived from geometry: {version}")

    mod_x = dist_x / (size - 7) if dist_x > 0 else unit
    mod_y = dist_y / (size - 7) if dist_y > 0 else unit
    origin_x = tl[0] - 3.5 * mod_x
    origin_y = tl[1] - 3.5 * mod_y

    def sample(r, c):
        x = int(round(origin_x + (c + 0.5) * mod_x))
        y = int(round(origin_y + (r + 0.5) * mod_y))
        return 1 if is_dark(x, y) else 0

    # --- format info: brute force against all 32 possibilities ---
    # raw bit positions for the first copy (matches encoder's set_format_info)
    def format_bits_from_grid():
        bits = [0] * 15
        for i in range(6):
            bits[i] = sample(i, 8)
        bits[6] = sample(7, 8)
        bits[7] = sample(8, 8)
        bits[8] = sample(8, 7)
        for i in range(9, 15):
            bits[i] = sample(8, 14 - i)
        val = 0
        for i in range(15):
            val |= bits[i] << i
        return val

    raw_fmt = format_bits_from_grid()
    best_fmt = None
    for ec_ind in range(4):
        for mask_id in range(8):
            fmt5 = (ec_ind << 3) | mask_id
            code = bch_format_info(fmt5)
            dist = bin(code ^ raw_fmt).count('1')
            if best_fmt is None or dist < best_fmt[0]:
                best_fmt = (dist, ec_ind, mask_id)
    if best_fmt is None or best_fmt[0] > 7:
        raise QRDecodeError("could not decode format info")
    _, ec_ind, mask_id = best_fmt
    level = EC_LEVEL_FROM_INDICATOR.get(ec_ind, 'L')
    if verbose:
        print(f"    detected version={version} size={size} level={level} mask={mask_id} "
              f"(format hamming dist={best_fmt[0]})")

    # --- read data bits following the encoder's zig-zag order ---
    skel = QRSkeleton(version).build()
    order = skel.data_order()
    bits = []
    for (r, c) in order:
        b = sample(r, c)
        if mask_condition(mask_id, r, c):
            b ^= 1
        bits.append(b)

    total_cw, ecc_cw, b1, dc1, b2, dc2 = QR_L[version]
    all_cw = bytearray()
    for i in range(0, min(len(bits), total_cw * 8), 8):
        byte = 0
        for j in range(8):
            byte = (byte << 1) | (bits[i + j] if i + j < len(bits) else 0)
        all_cw.append(byte)
    if len(all_cw) < total_cw:
        raise QRDecodeError("not enough codewords sampled from image")

    # --- de-interleave blocks (reverse of offline_qr_gen.py's interleave()) ---
    block_sizes = [dc1] * b1 + [dc2] * b2
    nblocks = b1 + b2
    max_data_len = max(block_sizes)
    data_blocks = [[] for _ in range(nblocks)]
    idx = 0
    for i in range(max_data_len):
        for bi in range(nblocks):
            if i < block_sizes[bi]:
                data_blocks[bi].append(all_cw[idx])
                idx += 1
    ecc_blocks = [[] for _ in range(nblocks)]
    for i in range(ecc_cw):
        for bi in range(nblocks):
            ecc_blocks[bi].append(all_cw[idx])
            idx += 1

    corrected_data = bytearray()
    for bi in range(nblocks):
        full_block = data_blocks[bi] + ecc_blocks[bi]
        try:
            corrected = rs_decode(full_block, ecc_cw)
        except RSDecodeError as e:
            raise QRDecodeError(f"block {bi}: {e}")
        corrected_data += bytes(corrected)

    # --- parse bit stream: mode + char count + byte data ---
    dbits = []
    for byte in corrected_data:
        for i in range(7, -1, -1):
            dbits.append((byte >> i) & 1)

    def read_bits(pos, n):
        v = 0
        for i in range(n):
            v = (v << 1) | (dbits[pos + i] if pos + i < len(dbits) else 0)
        return v

    mode = read_bits(0, 4)
    if mode != 0b0100:
        raise QRDecodeError(f"unsupported QR mode indicator {mode:#06b} (only byte mode supported)")
    ccbits = char_count_bits(version)
    length = read_bits(4, ccbits)
    data_start = 4 + ccbits
    out = bytearray()
    for i in range(length):
        out.append(read_bits(data_start + i * 8, 8))
    return bytes(out)


# =========================================================================
# PART i/N framing + CLI driver
# =========================================================================

def parse_part_payload(text):
    """Parse the 'PART<i>/<N>:<b64chunk>' framing used by offline_qr_gen.py."""
    if not text.startswith('PART'):
        return None
    rest = text[4:]
    if ':' not in rest:
        return None
    header, chunk = rest.split(':', 1)
    if '/' not in header:
        return None
    i_str, n_str = header.split('/', 1)
    try:
        i = int(i_str)
        n = int(n_str)
    except ValueError:
        return None
    return i, n, chunk


def main():
    ap = argparse.ArgumentParser(description="Decode PART i/N QR frames back into a file")
    ap.add_argument('--input-dir', '-i', required=True, help='directory of QR frame images (.ppm/.pgm/.pbm/.png)')
    ap.add_argument('--output', '-o', required=True, help='output file path to write the reconstructed data')
    ap.add_argument('--verbose', '-v', action='store_true')
    args = ap.parse_args()

    patterns = ('*.ppm', '*.pgm', '*.pbm', '*.png', '*.PPM', '*.PGM', '*.PBM', '*.PNG')
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(args.input_dir, pat)))
    files = sorted(set(files))
    if not files:
        print(f"No image files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    parts = {}
    expected_n = None
    for path in files:
        try:
            img = load_image(path)
            texts = [d.decode('ascii', errors='replace') for d in decode_qr_image_multi(img, verbose=args.verbose)]
        except Exception as e:
            print(f"[skip] {os.path.basename(path)}: {e}", file=sys.stderr)
            continue
        if args.verbose and len(texts) > 1:
            print(f"  {os.path.basename(path)}: decoded {len(texts)} QR code(s) from one image (grid mode)")
        for text in texts:
            parsed = parse_part_payload(text)
            if parsed is None:
                print(f"[skip] {os.path.basename(path)}: decoded but not a PART frame: {text[:60]!r}", file=sys.stderr)
                continue
            i, n, chunk = parsed
            if expected_n is None:
                expected_n = n
            elif n != expected_n:
                print(f"[warn] {os.path.basename(path)}: N={n} differs from expected {expected_n}, skipping", file=sys.stderr)
                continue
            if i in parts and parts[i] != chunk:
                print(f"[warn] {os.path.basename(path)}: duplicate part {i} with different content, keeping first", file=sys.stderr)
                continue
            if i not in parts:
                parts[i] = chunk
                if args.verbose:
                    print(f"  captured part {i}/{n} from {os.path.basename(path)}")

    if expected_n is None:
        print("No valid PART frames decoded from any image.", file=sys.stderr)
        sys.exit(1)

    missing = [i for i in range(1, expected_n + 1) if i not in parts]
    if missing:
        print(f"Missing parts: {missing} (have {len(parts)}/{expected_n})", file=sys.stderr)
        sys.exit(1)

    b64_text = ''.join(parts[i] for i in range(1, expected_n + 1))
    try:
        raw = base64.b64decode(b64_text)
    except Exception as e:
        print(f"base64 decode failed: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, 'wb') as f:
        f.write(raw)
    print(f"Wrote {len(raw)} bytes to {args.output} ({expected_n} part(s) reassembled)")


if __name__ == '__main__':
    main()
