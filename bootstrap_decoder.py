#!/usr/bin/env python3
"""Bootstrap grid-code decoder. Pure stdlib, no dependencies.

Type this file into a fresh .py file on the offline laptop, then run:
    python3 bootstrap_decoder.py <photos_dir> <output_file>
where <photos_dir> holds one .ppm photo per frame (see bootstrap_sender.html
for how the phone renders and how to save camera shots as .ppm).

GRID LAYOUT (GRID x GRID cells, same on sender and decoder):
  - Outer ring (row 0, row GRID-1, col 0, col GRID-1) is solid black.
    It is just a visual border to help you frame the photo; it is not read.
  - The inner (GRID-2) x (GRID-2) cells hold header + payload bits,
    filled row by row, left to right. Black cell = bit 1, white cell = bit 0.
  - Header = 48 bits (6 bytes), in this order:
        frame_index   (16 bits)
        total_frames  (16 bits)
        payload_len   (8 bits)   -- number of payload bytes in this frame
        checksum      (8 bits)   -- sum(payload bytes) mod 256
  - Payload = payload_len bytes, 8 bits each, most-significant bit first.
"""
import sys, os, glob

GRID = 32           # cells per side, including the black border ring
INNER = GRID - 2    # data cells per side
HEADER_BITS = 48


def read_ppm(path):
    with open(path, 'rb') as f:
        data = f.read()
    if not data.startswith(b'P6'):
        raise ValueError('not a P6 PPM file: ' + path)
    pos = 2
    vals = []
    while len(vals) < 3:
        while data[pos] in b' \t\r\n':
            pos += 1
        if data[pos:pos + 1] == b'#':
            while data[pos] not in b'\r\n':
                pos += 1
            continue
        start = pos
        while data[pos] not in b' \t\r\n':
            pos += 1
        vals.append(int(data[start:pos]))
    width, height, maxval = vals
    pos += 1  # single whitespace byte before the binary pixel data
    pixels = data[pos:pos + width * height * 3]
    return width, height, pixels


def cell_bit(pixels, width, height, row, col):
    """Average brightness of the middle half of one cell; dark => 1."""
    cw = width / GRID
    ch = height / GRID
    x0, x1 = int((col + 0.25) * cw), max(int((col + 0.75) * cw), 0)
    y0, y1 = int((row + 0.25) * ch), max(int((row + 0.75) * ch), 0)
    x1 = max(x1, x0 + 1)
    y1 = max(y1, y0 + 1)
    total = 0
    count = 0
    for y in range(y0, y1):
        row_off = y * width * 3
        for x in range(x0, x1):
            o = row_off + x * 3
            total += pixels[o] + pixels[o + 1] + pixels[o + 2]
            count += 1
    avg = total / (count * 3)
    return 1 if avg < 128 else 0


def bits_to_int(bits):
    value = 0
    for bit in bits:
        value = (value << 1) | bit
    return value


def decode_frame(path):
    width, height, pixels = read_ppm(path)
    bits = []
    for i in range(INNER * INNER):
        row = 1 + i // INNER
        col = 1 + i % INNER
        bits.append(cell_bit(pixels, width, height, row, col))

    header = bits[:HEADER_BITS]
    idx = bits_to_int(header[0:16])
    total = bits_to_int(header[16:32])
    plen = bits_to_int(header[32:40])
    checksum = bits_to_int(header[40:48])

    payload_bits = bits[HEADER_BITS:HEADER_BITS + plen * 8]
    if len(payload_bits) < plen * 8:
        return None  # frame too small / misread, can't trust it

    payload = bytearray()
    for b in range(plen):
        payload.append(bits_to_int(payload_bits[b * 8:(b + 1) * 8]))

    if sum(payload) % 256 != checksum:
        return None  # checksum mismatch, discard this frame

    return idx, total, bytes(payload)


def main():
    if len(sys.argv) != 3:
        print('usage: bootstrap_decoder.py <photos_dir> <output_file>')
        sys.exit(1)
    photos_dir, out_path = sys.argv[1], sys.argv[2]

    files = sorted(glob.glob(os.path.join(photos_dir, '*.ppm')))
    frames = {}
    total_expected = None
    skipped = 0

    for path in files:
        result = decode_frame(path)
        if result is None:
            skipped += 1
            print('skipped (bad checksum or unreadable):', path)
            continue
        idx, total, payload = result
        total_expected = total
        if idx in frames and frames[idx] != payload:
            print('warning: two different readings for frame', idx)
        frames[idx] = payload

    if total_expected is None:
        print('no valid frames decoded, nothing to write')
        sys.exit(1)

    missing = [i for i in range(total_expected) if i not in frames]
    if missing:
        print('missing frames, re-photograph these and re-run:', missing)
        sys.exit(1)

    with open(out_path, 'wb') as f:
        for i in range(total_expected):
            f.write(frames[i])

    print('OK: decoded %d frames (%d skipped) -> %s' % (len(frames), skipped, out_path))


if __name__ == '__main__':
    main()
