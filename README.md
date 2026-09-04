# QR Transfer — Airgapped File Transfer Kit

Move files to/from a laptop that has **no internet, no USB, and no LAN/WiFi
sharing**, using nothing but QR codes and a phone camera. Everything here is
pure Python 3 (standard library only — no `pip install` needed) plus a single
self-contained HTML file for the phone.

## Files

| File | Runs on | Needs internet/installs? | Purpose |
|---|---|---|---|
| `offline_qr_gen.py` | isolated laptop | No | Encodes a file into a cycling sequence of standard QR codes (terminal ASCII or PPM images) |
| `laptop_qr_decoder.py` | isolated laptop | No | Decodes photos of QR codes (single or grid) back into the original file |
| `phone_qr_scanner.html` | phone | No (once saved locally) | Scans QR codes cycling on the laptop screen and reassembles the file; also can display cycling QR codes to send a file *into* the laptop |
| `bootstrap_decoder.py` | isolated laptop | No | Tiny (~150 line) fallback decoder for the *very first* transfer, small enough to type by hand if you have no way to get the two files above onto the laptop at all |
| `bootstrap_sender.html` | phone | No (once saved locally) | Companion sender for `bootstrap_decoder.py`'s simplified grid-code format |

## Prerequisites

- The isolated laptop already has Python 3 installed.
- You have a phone with a camera and (at some point, before you need it
  offline) normal internet access to save `phone_qr_scanner.html` locally.

**Get `phone_qr_scanner.html` onto your phone now**, while you have internet —
email it to yourself, AirDrop it, etc. Once saved, it works fully offline.

## One-time bootstrap (only if you have literally no way to copy files onto the laptop)

If `offline_qr_gen.py` / `laptop_qr_decoder.py` aren't already on the isolated
laptop, and there is truly no USB/Bluetooth/SD/CD/printer-scanner channel
available, use the bootstrap:

1. On the laptop, hand-type `bootstrap_decoder.py` into a new file exactly as
   written, and save it. It's short and dependency-free by design.
2. Get `bootstrap_sender.html` onto your phone (email/AirDrop, same as above).
3. On the phone, open `bootstrap_sender.html`, choose `offline_qr_gen.py` as
   the file to send, and start cycling.
4. On the laptop, photograph every cycling frame with its built-in camera app,
   saving each shot as a `.ppm` file into one folder (the decoder only reads
   `.ppm` — you'll need some existing way to get camera output into that
   format; this is the main practical limitation of the bootstrap path).
5. Run:
   ```bash
   python3 bootstrap_decoder.py <photos_folder> offline_qr_gen.py
   ```
   If it reports missing/bad frames, re-photograph just those and re-run.
6. Repeat steps 3-5 for `laptop_qr_decoder.py`.

Note: for a ~28KB file this scheme needs on the order of 250+ photos, since
each frame only carries about 106 bytes. Use this only as a last resort —
prefer any other physical channel (Bluetooth, SD card, CD/DVD, printer +
scanner/OCR) if one is actually available, since it will be far faster.

## Normal use, once the two main scripts are on the laptop

### Get files OUT of the isolated laptop (laptop → phone)

```bash
# on the laptop
zip -r mycode.zip ./my-project
python3 offline_qr_gen.py mycode.zip
```

On the phone, open `phone_qr_scanner.html` → **Scan** mode → **Start Camera** →
point at the laptop screen. Progress shows as "X/N parts captured"; once
complete it auto-downloads the reconstructed file.

### Get files INTO the isolated laptop (phone → laptop)

On the phone, open `phone_qr_scanner.html` → **Send** mode → pick a **grid
size** (e.g. `2x2` to cut the number of photos roughly 4x) → choose the file →
**Start Cycling**.

On the laptop, photograph each displayed frame with its camera app, saving
into one folder, then:

```bash
python3 laptop_qr_decoder.py -i <photos_folder> -o received_file
```

Missing or unreadable parts are reported by number — re-photograph just those
and re-run the same command.

### Command reference

```bash
# Generate QR frames in the terminal (default, single QR per frame)
python3 offline_qr_gen.py <file>

# Generate QR frames as PPM images instead of terminal display
python3 offline_qr_gen.py <file> --ppm-dir out/

# Pack multiple QR codes per frame/photo (fewer photos needed)
python3 offline_qr_gen.py <file> --grid 2x2 --ppm-dir out/

# Custom cycle delay (seconds) for terminal display
python3 offline_qr_gen.py <file> --delay 2.0

# Single pass instead of looping forever (useful for scripting/testing)
python3 offline_qr_gen.py <file> --once

# Decode a folder of QR photos (single or grid images, auto-detected)
python3 laptop_qr_decoder.py -i <photos_folder> -o <output_file>
python3 laptop_qr_decoder.py -i <photos_folder> -o <output_file> -v   # verbose
```

## How it works

Files are base64-encoded and split into chunks, each wrapped as
`PART<i>/<N>:<chunk>` and rendered as one standard QR code (error-correction
level L). The receiving side collects parts by index (order-independent,
duplicate-tolerant), and once all `N` parts are present, concatenates and
base64-decodes them back into the original bytes. `laptop_qr_decoder.py`
implements QR image decoding entirely from scratch (finder-pattern detection,
Reed-Solomon error correction, PNG/PPM parsing) — no external libraries.

## Known limitations

- Photo-based decoding needs clean, well-lit, head-on, non-blurry shots — it's
  a best-effort tool, not a commercial-grade scanner.
- Higher `--grid` sizes trade capture speed for reliability (smaller QR codes
  per tile); drop to a smaller grid if scans keep failing.
- The bootstrap path only reads `.ppm` images and needs several hundred photos
  for a real-sized file — treat it strictly as a last-resort, one-time
  mechanism, not a regular workflow.
