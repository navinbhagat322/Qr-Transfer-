#!/usr/bin/env python3
"""Self-contained QR Code encoder (byte mode, EC level L) + CLI transfer tool.

Implements QR encoding from scratch: GF(256) math, Reed-Solomon ECC,
version selection, format/version info BCH, byte-mode encoding, all 8
mask patterns with penalty scoring, and module placement.

CLI splits a file's base64 into chunks, each wrapped as "PART<i>/<N>:<data>"
byte-mode QR symbol, then renders as terminal ASCII (looping) or PPM files.
"""
import sys
import os
import io
import time
import base64
import zipfile
import struct
import argparse

# ----------------------------- GF(256) ----------------------------------

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

# ------------------------- Reed-Solomon ----------------------------------

def rs_generator_poly(n_ecc):
    """Generator polynomial coefficients (highest degree first), monic."""
    g = [1]
    for i in range(n_ecc):
        g.append(0)
        for j in range(len(g) - 1, 0, -1):
            g[j] ^= gf_mul(g[j - 1], _EXP[i])
    return g

def rs_encode(data, n_ecc):
    gen = rs_generator_poly(n_ecc)
    res = list(data) + [0] * n_ecc
    for i in range(len(data)):
        coef = res[i]
        if coef == 0:
            continue
        for j in range(len(gen)):
            res[i + j] ^= gf_mul(gen[j], coef)
    return res[len(data):]

# ------------------------- QR tables (EC level L) -------------------------
# Per-version: (total_codewords, ecc_codewords_per_block, num_blocks_group1,
#               data_codewords_per_block_group1, num_blocks_group2, data_codewords_per_block_group2)
# Level L data from the QR spec, versions 1..40.
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

# Character count indicator bit lengths for byte mode by version range
def char_count_bits(version):
    if version <= 9:
        return 8
    elif version <= 26:
        return 16
    else:
        return 16

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

REMAINDER_BITS = [0,0,7,7,7,7,7,0,0,0,0,0,0,0,3,3,3,3,3,3,3,4,4,4,4,4,4,4,3,3,3,3,3,3,3,0,0,0,0,0,0]

# ---------------------- BCH format/version info ---------------------------

def bch_format_info(fmt5):
    """fmt5: 5 bits (2 EC level bits + 3 mask bits). Return 15-bit format info."""
    g = 0b10100110111  # 0x537, degree 10
    data = fmt5 << 10
    val = data
    while val.bit_length() > 10:
        val ^= g << (val.bit_length() - 11)
    fmt = (data | val) ^ 0b101010000010010  # mask
    return fmt

def bch_version_info(version):
    """Return 18-bit version info for versions 7-40."""
    g = 0b1111100100101  # 0x1F25, degree 12
    data = version << 12
    val = data
    while val.bit_length() > 12:
        val ^= g << (val.bit_length() - 13)
    return data | val

EC_LEVEL_INDICATOR = {'L': 0b01, 'M': 0b00, 'Q': 0b11, 'H': 0b10}

# ------------------------------ Bit buffer ---------------------------------

class BitBuffer:
    def __init__(self):
        self.bits = []
    def put(self, val, length):
        for i in range(length - 1, -1, -1):
            self.bits.append((val >> i) & 1)
    def __len__(self):
        return len(self.bits)
    def to_bytes(self):
        out = bytearray()
        for i in range(0, len(self.bits), 8):
            byte = 0
            for j in range(8):
                bit = self.bits[i + j] if i + j < len(self.bits) else 0
                byte = (byte << 1) | bit
            out.append(byte)
        return bytes(out)

# ------------------------------ Version selection --------------------------

def choose_version(data_len, level='L'):
    """Choose smallest version whose data capacity (bytes) fits byte-mode encoded data_len bytes."""
    for version in range(1, 41):
        total_cw, ecc_cw, b1, dc1, b2, dc2 = QR_L[version]
        data_capacity_cw = b1 * dc1 + b2 * dc2
        ccbits = char_count_bits(version)
        header_bits = 4 + ccbits
        needed_bits = header_bits + data_len * 8
        # plus terminator up to 4 bits, but capacity check is fine w/ approx
        needed_bytes = (needed_bits + 4 + 7) // 8
        if needed_bytes <= data_capacity_cw:
            return version
    raise ValueError("data too large for QR (even version 40, level L)")

def max_bytes_for_version(version):
    total_cw, ecc_cw, b1, dc1, b2, dc2 = QR_L[version]
    data_capacity_cw = b1 * dc1 + b2 * dc2
    ccbits = char_count_bits(version)
    header_bits = 4 + ccbits
    avail_bits = data_capacity_cw * 8 - header_bits
    return avail_bits // 8

# ------------------------------ Encoding ------------------------------------

def encode_data_codewords(data_bytes, version):
    total_cw, ecc_cw, b1, dc1, b2, dc2 = QR_L[version]
    data_capacity_cw = b1 * dc1 + b2 * dc2
    ccbits = char_count_bits(version)

    bb = BitBuffer()
    bb.put(0b0100, 4)  # byte mode
    bb.put(len(data_bytes), ccbits)
    for byte in data_bytes:
        bb.put(byte, 8)

    # terminator
    remaining = data_capacity_cw * 8 - len(bb)
    bb.put(0, min(4, max(0, remaining)))

    # pad to byte boundary
    while len(bb) % 8 != 0:
        bb.bits.append(0)

    codewords = bytearray(bb.to_bytes())

    # pad bytes
    pad_bytes = [0xEC, 0x11]
    i = 0
    while len(codewords) < data_capacity_cw:
        codewords.append(pad_bytes[i % 2])
        i += 1

    return bytes(codewords[:data_capacity_cw])

def build_blocks(data_codewords, version):
    total_cw, ecc_cw, b1, dc1, b2, dc2 = QR_L[version]
    blocks = []
    idx = 0
    for _ in range(b1):
        blocks.append(data_codewords[idx:idx + dc1])
        idx += dc1
    for _ in range(b2):
        blocks.append(data_codewords[idx:idx + dc2])
        idx += dc2
    ecc_blocks = [rs_encode(block, ecc_cw) for block in blocks]
    return blocks, ecc_blocks

def interleave(blocks, ecc_blocks):
    max_data_len = max(len(b) for b in blocks)
    result = []
    for i in range(max_data_len):
        for b in blocks:
            if i < len(b):
                result.append(b[i])
    ecc_len = len(ecc_blocks[0])
    for i in range(ecc_len):
        for eb in ecc_blocks:
            result.append(eb[i])
    return bytes(result)

# ------------------------------ Matrix building ------------------------------

class QRMatrix:
    def __init__(self, version):
        self.version = version
        self.size = version * 4 + 17
        self.modules = [[None] * self.size for _ in range(self.size)]  # None=unset, 0/1
        self.is_function = [[False] * self.size for _ in range(self.size)]

    def set(self, r, c, val, function=True):
        if 0 <= r < self.size and 0 <= c < self.size:
            self.modules[r][c] = val
            if function:
                self.is_function[r][c] = True

    def place_finder(self, r, c):
        for dr in range(-1, 8):
            for dc in range(-1, 8):
                rr, cc = r + dr, c + dc
                if not (0 <= rr < self.size and 0 <= cc < self.size):
                    continue
                if 0 <= dr <= 6 and 0 <= dc <= 6:
                    is_border = dr in (0, 6) or dc in (0, 6)
                    is_inner = 2 <= dr <= 4 and 2 <= dc <= 4
                    val = 1 if (is_border or is_inner) else 0
                    self.set(rr, cc, val)
                else:
                    self.set(rr, cc, 0)  # separator

    def place_finders(self):
        self.place_finder(0, 0)
        self.place_finder(0, self.size - 7)
        self.place_finder(self.size - 7, 0)

    def place_alignment(self):
        coords = ALIGNMENT_COORDS[self.version]
        for r in coords:
            for c in coords:
                # skip if overlapping finder patterns
                if (r <= 8 and c <= 8) or (r <= 8 and c >= self.size - 9) or (r >= self.size - 9 and c <= 8):
                    continue
                for dr in range(-2, 3):
                    for dc in range(-2, 3):
                        d = max(abs(dr), abs(dc))
                        val = 0 if d == 1 else 1
                        self.set(r + dr, c + dc, val)

    def place_timing(self):
        for i in range(8, self.size - 8):
            val = 1 if i % 2 == 0 else 0
            if self.modules[6][i] is None:
                self.set(6, i, val)
            if self.modules[i][6] is None:
                self.set(i, 6, val)

    def place_dark_module(self):
        self.set(self.size - 8, 8, 1)

    def reserve_format_areas(self):
        for i in range(9):
            if self.modules[8][i] is None:
                self.set(8, i, 0)
            if self.modules[i][8] is None:
                self.set(i, 8, 0)
        for i in range(self.size - 8, self.size):
            self.set(8, i, 0)
        for i in range(self.size - 7, self.size):
            self.set(i, 8, 0)

    def reserve_version_areas(self):
        if self.version < 7:
            return
        for r in range(6):
            for c in range(self.size - 11, self.size - 8):
                self.set(r, c, 0)
        for c in range(6):
            for r in range(self.size - 11, self.size - 8):
                self.set(r, c, 0)

    def place_data(self, data_bits):
        bit_idx = 0
        n = len(data_bits)
        col = self.size - 1
        upward = True
        while col > 0:
            if col == 6:
                col -= 1
            for i in range(self.size):
                row = (self.size - 1 - i) if upward else i
                for c in (col, col - 1):
                    if not self.is_function[row][c]:
                        bit = data_bits[bit_idx] if bit_idx < n else 0
                        self.modules[row][c] = bit
                        bit_idx += 1
            upward = not upward
            col -= 2

    def apply_mask(self, mask_id, only_data=True):
        for r in range(self.size):
            for c in range(self.size):
                if only_data and self.is_function[r][c]:
                    continue
                if mask_condition(mask_id, r, c):
                    self.modules[r][c] ^= 1

    def set_format_info(self, ec_indicator, mask_id):
        fmt5 = (ec_indicator << 3) | mask_id
        bits = bch_format_info(fmt5)  # 15-bit value, bit0 = LSB

        def gb(i):
            return (bits >> i) & 1

        # first copy
        for i in range(6):
            self.set(i, 8, gb(i))
        self.set(7, 8, gb(6))
        self.set(8, 8, gb(7))
        self.set(8, 7, gb(8))
        for i in range(9, 15):
            self.set(8, 14 - i, gb(i))
        # second copy
        for i in range(8):
            self.set(8, self.size - 1 - i, gb(i))
        for i in range(8, 15):
            self.set(self.size - 15 + i, 8, gb(i))

    def set_version_info(self):
        if self.version < 7:
            return
        bits = bch_version_info(self.version)  # 18-bit value, bit0 = LSB
        for i in range(18):
            b = (bits >> i) & 1
            a = self.size - 11 + i % 3
            row = i // 3
            self.set(row, a, b)
            self.set(a, row, b)

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

# --------------------------- Penalty scoring ---------------------------------

def penalty_score(mat):
    size = mat.size
    m = mat.modules
    score = 0
    # rule 1: consecutive same-color modules in row/col
    for r in range(size):
        run = 1
        for c in range(1, size):
            if m[r][c] == m[r][c - 1]:
                run += 1
            else:
                if run >= 5:
                    score += 3 + (run - 5)
                run = 1
        if run >= 5:
            score += 3 + (run - 5)
    for c in range(size):
        run = 1
        for r in range(1, size):
            if m[r][c] == m[r - 1][c]:
                run += 1
            else:
                if run >= 5:
                    score += 3 + (run - 5)
                run = 1
        if run >= 5:
            score += 3 + (run - 5)
    # rule 2: 2x2 blocks same color
    for r in range(size - 1):
        for c in range(size - 1):
            v = m[r][c]
            if v == m[r][c + 1] == m[r + 1][c] == m[r + 1][c + 1]:
                score += 3
    # rule 3: finder-like patterns 1:1:3:1:1 with 4 white either side
    pattern = [1, 0, 1, 1, 1, 0, 1]
    pattern_ext_a = [0, 0, 0, 0] + pattern
    pattern_ext_b = pattern + [0, 0, 0, 0]
    for r in range(size):
        row = m[r]
        for c in range(size - 6):
            seg = row[c:c + 7]
            if seg == pattern:
                if c >= 4 and row[c - 4:c + 7] == pattern_ext_a:
                    score += 40
                elif c + 11 <= size and row[c:c + 11] == pattern_ext_b:
                    score += 40
    for c in range(size):
        col = [m[r][c] for r in range(size)]
        for r in range(size - 6):
            seg = col[r:r + 7]
            if seg == pattern:
                if r >= 4 and col[r - 4:r + 7] == pattern_ext_a:
                    score += 40
                elif r + 11 <= size and col[r:r + 11] == pattern_ext_b:
                    score += 40
    # rule 4: proportion of dark modules
    dark = sum(sum(row) for row in m)
    total = size * size
    percent = dark * 100 / total
    prev_mult = int(percent // 5) * 5
    next_mult = prev_mult + 5
    a = abs(prev_mult - 50) // 5
    b = abs(next_mult - 50) // 5
    score += min(a, b) * 10
    return score

# --------------------------- Top-level QR build ------------------------------

def build_qr(data_bytes, version=None, level='L'):
    if version is None:
        version = choose_version(len(data_bytes), level)
    data_cw = encode_data_codewords(data_bytes, version)
    blocks, ecc_blocks = build_blocks(data_cw, version)
    all_cw = interleave(blocks, ecc_blocks)
    bits = []
    for byte in all_cw:
        for i in range(7, -1, -1):
            bits.append((byte >> i) & 1)
    remainder = REMAINDER_BITS[version - 1]
    bits.extend([0] * remainder)

    best = None
    ec_ind = EC_LEVEL_INDICATOR[level]
    for mask_id in range(8):
        mat = QRMatrix(version)
        mat.place_finders()
        mat.place_alignment()
        mat.place_timing()
        mat.place_dark_module()
        mat.reserve_format_areas()
        mat.reserve_version_areas()
        mat.place_data(bits)
        mat.apply_mask(mask_id)
        mat.set_format_info(ec_ind, mask_id)
        mat.set_version_info()
        score = penalty_score(mat)
        if best is None or score < best[0]:
            best = (score, mat)
    return best[1]

# ------------------------------ Rendering ------------------------------------

def render_ascii(mat, quiet=4):
    size = mat.size
    lines = []
    total = size + 2 * quiet
    blank_row = "  " * total
    for _ in range(quiet):
        lines.append(blank_row)
    for r in range(size):
        row_chars = []
        for _ in range(quiet):
            row_chars.append("  ")
        for c in range(size):
            row_chars.append("██" if mat.modules[r][c] else "  ")
        for _ in range(quiet):
            row_chars.append("  ")
        lines.append("".join(row_chars))
    for _ in range(quiet):
        lines.append(blank_row)
    return "\n".join(lines)

def parse_grid_spec(s):
    """Parse a 'ROWSxCOLS' grid spec, e.g. '2x2' or '3x2'."""
    s = s.lower().strip()
    if 'x' not in s:
        raise ValueError("grid spec must look like ROWSxCOLS, e.g. 2x2")
    r_str, c_str = s.split('x', 1)
    rows, cols = int(r_str), int(c_str)
    if rows < 1 or cols < 1:
        raise ValueError("grid rows/cols must be >= 1")
    return rows, cols


def chunk_into_grids(items, rows, cols):
    """Split a flat list into groups of size rows*cols (last group padded with None)."""
    per_grid = rows * cols
    groups = []
    for i in range(0, len(items), per_grid):
        group = items[i:i + per_grid]
        while len(group) < per_grid:
            group.append(None)
        groups.append(group)
    return groups


def render_grid_ascii(cell_mats, rows, cols, quiet=2, gap=2):
    """Render a rows x cols grid of QRMatrix objects (None = empty cell) as one
    combined ASCII screen. Each cell keeps its own quiet zone; cells are
    separated by `gap` blank module-columns/rows."""
    cell_texts = []
    max_h = 0
    max_w = 0
    for mat in cell_mats:
        if mat is None:
            cell_texts.append(None)
            continue
        lines = render_ascii(mat, quiet=quiet).split("\n")
        cell_texts.append(lines)
        max_h = max(max_h, len(lines))
        max_w = max(max_w, len(lines[0]) if lines else 0)

    gap_str = "  " * gap
    out_rows = []
    for r in range(rows):
        # build max_h lines for this row of cells
        row_lines = [""] * max_h
        for c in range(cols):
            idx = r * cols + c
            lines = cell_texts[idx] if idx < len(cell_texts) else None
            for li in range(max_h):
                if lines is not None and li < len(lines):
                    piece = lines[li].ljust(max_w)
                else:
                    piece = " " * max_w
                sep = gap_str if c > 0 else ""
                row_lines[li] += sep + piece
        out_rows.extend(row_lines)
        if r < rows - 1:
            out_rows.append("")
    return "\n".join(out_rows)


def render_grid_ppm(cell_mats, rows, cols, path, scale=8, quiet=4, gap=2):
    """Render a rows x cols grid of QRMatrix objects (None = empty cell) into
    a single combined PPM image. Each cell is allocated a uniform-size box
    (based on the largest QR version present) with its own quiet zone,
    separated by `gap` blank modules of white space."""
    present = [m for m in cell_mats if m is not None]
    if not present:
        raise ValueError("grid has no QR codes to render")
    max_size = max(m.size for m in present)
    cell_dim = (max_size + 2 * quiet) * scale
    gap_px = gap * scale
    img_w = cols * cell_dim + (cols - 1) * gap_px
    img_h = rows * cell_dim + (rows - 1) * gap_px

    buf = bytearray(b'\xff' * (img_w * img_h * 3))
    white_scale_row = b'\xff' * (scale * 3)
    black_scale_row = b'\x00' * (scale * 3)

    for idx, mat in enumerate(cell_mats):
        if mat is None:
            continue
        r_idx, c_idx = divmod(idx, cols)
        ox = c_idx * (cell_dim + gap_px)
        oy = r_idx * (cell_dim + gap_px)
        size = mat.size
        for rr in range(size):
            py0 = oy + (quiet + rr) * scale
            for cc in range(size):
                px0 = ox + (quiet + cc) * scale
                block = black_scale_row if mat.modules[rr][cc] else white_scale_row
                for py in range(py0, py0 + scale):
                    rowstart = (py * img_w + px0) * 3
                    buf[rowstart:rowstart + scale * 3] = block

    header = f"P6\n{img_w} {img_h}\n255\n".encode('ascii')
    with open(path, 'wb') as f:
        f.write(header)
        f.write(bytes(buf))


def render_ppm(mat, path, scale=8, quiet=4):
    size = mat.size
    total = size + 2 * quiet
    img_size = total * scale
    header = f"P6\n{img_size} {img_size}\n255\n".encode('ascii')
    row_bytes = bytearray()
    black = bytes([0, 0, 0])
    white = bytes([255, 255, 255])
    with open(path, 'wb') as f:
        f.write(header)
        for r in range(img_size):
            mod_row = r // scale - quiet
            line = bytearray()
            for c in range(img_size):
                mod_col = c // scale - quiet
                if 0 <= mod_row < size and 0 <= mod_col < size and mat.modules[mod_row][mod_col]:
                    line += black
                else:
                    line += white
            f.write(bytes(line))

# ------------------------------ CLI logic -------------------------------------

def clear_screen():
    if os.name == 'nt':
        os.system('cls')
    else:
        sys.stdout.write('\x1b[2J\x1b[H')
        sys.stdout.flush()

def read_input_bytes(path):
    if os.path.isdir(path):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
            for root, dirs, files in os.walk(path):
                for name in files:
                    full = os.path.join(root, name)
                    arcname = os.path.relpath(full, path)
                    zf.write(full, arcname)
        return buf.getvalue()
    else:
        with open(path, 'rb') as f:
            return f.read()

def build_chunks(b64_text, target_version, level='L'):
    # reserve header overhead: "PART<i>/<N>:" -- estimate max header length
    # Try increasing N guess iteratively to fit header size correctly.
    max_bytes = max_bytes_for_version(target_version)
    # binary-search-ish: assume N up to 4 digits initially, refine
    n_guess = 1
    while True:
        header_len_guess = len(f"PART{n_guess}/{n_guess}:")
        chunk_data_len = max_bytes - header_len_guess
        if chunk_data_len <= 0:
            raise ValueError("version too small to fit header")
        n_chunks = max(1, -(-len(b64_text) // chunk_data_len))
        header_len_actual = len(f"PART{n_chunks}/{n_chunks}:")
        chunk_data_len_actual = max_bytes - header_len_actual
        n_chunks_actual = max(1, -(-len(b64_text) // chunk_data_len_actual))
        if n_chunks_actual == n_chunks:
            n_chunks = n_chunks_actual
            chunk_data_len = chunk_data_len_actual
            break
        n_guess = n_chunks_actual
    chunks = []
    idx = 0
    for i in range(1, n_chunks + 1):
        piece = b64_text[idx: idx + chunk_data_len]
        idx += chunk_data_len
        chunks.append(f"PART{i}/{n_chunks}:{piece}")
    return chunks

def main():
    ap = argparse.ArgumentParser(description="Offline file-transfer QR generator")
    ap.add_argument('input_file')
    ap.add_argument('--version', type=int, default=20, help='target QR version for chunk sizing (default 20)')
    ap.add_argument('--level', default='L', choices=['L'], help='EC level (only L implemented)')
    ap.add_argument('--delay', type=float, default=1.5, help='seconds between frames in loop mode')
    ap.add_argument('--once', '--no-loop', dest='once', action='store_true', help='render once and exit')
    ap.add_argument('--ppm-dir', default=None, help='write PPM frames to this directory instead of / in addition to looping')
    ap.add_argument('--scale', type=int, default=8, help='pixels per module for PPM output')
    ap.add_argument('--grid', default=None, help="render ROWSxCOLS QR codes together per screen/image, e.g. '2x2' or '3x2' (default: single QR per frame)")
    args = ap.parse_args()

    raw = read_input_bytes(args.input_file)
    b64_text = base64.b64encode(raw).decode('ascii')
    chunks = build_chunks(b64_text, args.version, args.level)
    n = len(chunks)

    frames = []
    for chunk_str in chunks:
        data_bytes = chunk_str.encode('ascii')
        version = choose_version(len(data_bytes), args.level)
        mat = build_qr(data_bytes, version=version, level=args.level)
        frames.append((mat, chunk_str))

    grid = None
    if args.grid:
        grid = parse_grid_spec(args.grid)

    if grid:
        rows, cols = grid
        mats_only = [mat for mat, _ in frames]
        grid_groups = chunk_into_grids(mats_only, rows, cols)
        n_grids = len(grid_groups)

        if args.ppm_dir:
            os.makedirs(args.ppm_dir, exist_ok=True)
            for gi, group in enumerate(grid_groups, 1):
                path = os.path.join(args.ppm_dir, f"grid_{gi:03d}_of_{n_grids:03d}.ppm")
                render_grid_ppm(group, rows, cols, path, scale=args.scale)
            print(f"Wrote {n_grids} grid PPM image(s) ({rows}x{cols} = up to {rows*cols} QR codes each, "
                  f"{n} total parts) to {args.ppm_dir}")

        if args.once:
            for gi, group in enumerate(grid_groups, 1):
                present = sum(1 for m in group if m is not None)
                print(f"Grid {gi}/{n_grids} ({present} QR code(s), {rows}x{cols} layout)")
                print(render_grid_ascii(group, rows, cols))
            return

        if not args.ppm_dir:
            try:
                while True:
                    for gi, group in enumerate(grid_groups, 1):
                        clear_screen()
                        present = sum(1 for m in group if m is not None)
                        print(f"Grid {gi}/{n_grids} ({present} QR code(s), {rows}x{cols} layout, {n} parts total)")
                        print(render_grid_ascii(group, rows, cols))
                        time.sleep(args.delay)
            except KeyboardInterrupt:
                print("\nStopped.")
        return

    if args.ppm_dir:
        os.makedirs(args.ppm_dir, exist_ok=True)
        for i, (mat, chunk_str) in enumerate(frames, 1):
            path = os.path.join(args.ppm_dir, f"part_{i:03d}_of_{n:03d}.ppm")
            render_ppm(mat, path, scale=args.scale)
        print(f"Wrote {n} PPM frame(s) to {args.ppm_dir}")

    if args.once:
        for i, (mat, chunk_str) in enumerate(frames, 1):
            print(f"Frame {i}/{n} (version {mat.version}, {len(chunk_str)} chars)")
            print(render_ascii(mat))
        return

    if not args.ppm_dir:
        try:
            while True:
                for i, (mat, chunk_str) in enumerate(frames, 1):
                    clear_screen()
                    print(f"Frame {i}/{n} (version {mat.version})")
                    print(render_ascii(mat))
                    time.sleep(args.delay)
        except KeyboardInterrupt:
            print("\nStopped.")

if __name__ == '__main__':
    main()
