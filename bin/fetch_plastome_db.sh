#!/bin/bash
#
# fetch_plastome_db.sh
#
# Download and install the Fast-Plast reference plastome database from Zenodo.
# The database is too large to ship in the git repository, so it is hosted as a
# versioned, checksummed archive and fetched on demand. The archive contains:
#     GenBank_Plastomes                  (FASTA, bare-accession headers)
#     GenBank_Plastomes.metadata.tsv     (accession -> taxonomy)
#
# After uploading a new database version to Zenodo, update DB_URL and DB_SHA256
# below (and DB_VERSION for the log). The checksum is what makes an install
# reproducible: everyone who runs this gets byte-identical reference data.
#
# Usage:
#     bin/fetch_plastome_db.sh [target_dir]
# target_dir defaults to the directory this script lives in (Fast-Plast/bin).

set -euo pipefail

# ---- configure per release -------------------------------------------------
DB_VERSION="v2.0"
# Zenodo direct-download URL. Format:
#   https://zenodo.org/records/<RECORD_ID>/files/<FILENAME>?download=1
DB_URL="https://zenodo.org/records/RECORD_ID/files/GenBank_Plastomes_v2.tar.gz?download=1"
DB_SHA256="PUT_REAL_SHA256_HERE"
ARCHIVE_NAME="GenBank_Plastomes_v2.tar.gz"
# Files expected inside the archive (used to verify a good extraction):
EXPECTED_FASTA="GenBank_Plastomes"
EXPECTED_META="GenBank_Plastomes.metadata.tsv"
# ---------------------------------------------------------------------------

TARGET_DIR="${1:-$(cd "$(dirname "$0")" && pwd)}"
mkdir -p "$TARGET_DIR"

echo "Fast-Plast reference database $DB_VERSION"
echo "  target: $TARGET_DIR"

# already installed and intact? then we're done.
if [[ -s "$TARGET_DIR/$EXPECTED_FASTA" && -s "$TARGET_DIR/$EXPECTED_META" ]]; then
    echo "  database already present; nothing to do."
    echo "  (delete $TARGET_DIR/$EXPECTED_FASTA to force a re-download)"
    exit 0
fi

# pick a downloader
if command -v curl >/dev/null 2>&1; then
    DL() { curl -L --fail --retry 3 -o "$1" "$2"; }
elif command -v wget >/dev/null 2>&1; then
    DL() { wget -O "$1" "$2"; }
else
    echo "ERROR: need either curl or wget to download the database." >&2
    exit 1
fi

# pick a sha256 tool (Linux: sha256sum, macOS: shasum -a 256)
if command -v sha256sum >/dev/null 2>&1; then
    SHA() { sha256sum "$1" | awk '{print $1}'; }
elif command -v shasum >/dev/null 2>&1; then
    SHA() { shasum -a 256 "$1" | awk '{print $1}'; }
else
    echo "ERROR: need sha256sum or shasum to verify the download." >&2
    exit 1
fi

if [[ "$DB_URL" == *RECORD_ID* || "$DB_SHA256" == "PUT_REAL_SHA256_HERE" ]]; then
    echo "ERROR: this script has not been configured for a release yet." >&2
    echo "       Edit DB_URL and DB_SHA256 in $0 after uploading to Zenodo." >&2
    exit 1
fi

ARCHIVE_PATH="$TARGET_DIR/$ARCHIVE_NAME"
echo "  downloading $ARCHIVE_NAME ..."
DL "$ARCHIVE_PATH" "$DB_URL"

echo "  verifying checksum ..."
GOT="$(SHA "$ARCHIVE_PATH")"
if [[ "$GOT" != "$DB_SHA256" ]]; then
    echo "ERROR: checksum mismatch -- download is corrupt or the wrong file." >&2
    echo "  expected: $DB_SHA256" >&2
    echo "  got:      $GOT" >&2
    echo "  (leaving $ARCHIVE_PATH in place for inspection)" >&2
    exit 1
fi
echo "  checksum OK."

echo "  extracting ..."
tar -xzf "$ARCHIVE_PATH" -C "$TARGET_DIR"

if [[ ! -s "$TARGET_DIR/$EXPECTED_FASTA" || ! -s "$TARGET_DIR/$EXPECTED_META" ]]; then
    echo "ERROR: archive extracted but expected files are missing." >&2
    echo "       wanted: $EXPECTED_FASTA and $EXPECTED_META" >&2
    exit 1
fi

rm -f "$ARCHIVE_PATH"
echo "  installed:"
echo "    $TARGET_DIR/$EXPECTED_FASTA"
echo "    $TARGET_DIR/$EXPECTED_META"
echo "Done."
