#!/usr/bin/env bash
# Deploy the pinned official Linux x86_64 static binary within the project.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
DEST="$ROOT/bin/HitSV/2.0.4"
URL=https://github.com/hitbc/HitSV/releases/download/2.0.4/HitSV_2.0.4
SHA=a22d3a87d4f2cdcf8d7d06558583e3600a2539eb427f6fd35c32dc8def7a4c6e
COMMIT=98f344afbe0c453c7f73091e58ed0ba00711337f
SOURCE_SHA=6fd053a22783a400e319ddbd0b613f914ee90b0786ebe4bce66e2b27c920e097
SOURCE_URL="https://codeload.github.com/hitbc/HitSV/tar.gz/$COMMIT"
SOURCE_ARCHIVE="$DEST/HitSV-2.0.4-source.tar.gz"
PROFILE_SHA=b91da787e60a683369518c016069e15f000a2daf6827a6dac1354bcccef9b824
mkdir -p "$DEST"
if ! printf '%s  %s\n' "$SHA" "$DEST/HitSV" | sha256sum -c --status; then
    curl --fail --location --retry 2 "$URL" -o "$DEST/HitSV.download"
    printf '%s  %s\n' "$SHA" "$DEST/HitSV.download" | sha256sum -c
    chmod 755 "$DEST/HitSV.download"
    mv "$DEST/HitSV.download" "$DEST/HitSV"
fi
if ! printf '%s  %s\n' "$SOURCE_SHA" "$SOURCE_ARCHIVE" | sha256sum -c --status; then
    curl --fail --location --retry 2 "$SOURCE_URL" -o "$SOURCE_ARCHIVE.download"
    printf '%s  %s\n' "$SOURCE_SHA" "$SOURCE_ARCHIVE.download" | sha256sum -c
    mv "$SOURCE_ARCHIVE.download" "$SOURCE_ARCHIVE"
fi
# Official k-mer error profile used as the Hybrid -k default.
tar -xOf "$SOURCE_ARCHIVE" "HitSV-$COMMIT/demo/9mer_profile" > "$DEST/9mer_profile"
printf '%s  %s\n' "$PROFILE_SHA" "$DEST/9mer_profile" | sha256sum -c --status
printf '%s  HitSV\n%s  HitSV-2.0.4-source.tar.gz\n%s  9mer_profile\n' "$SHA" "$SOURCE_SHA" "$PROFILE_SHA" > "$DEST/SHA256SUMS"
printf 'Version: 2.0.4\nBinary: %s\nRelease: https://github.com/hitbc/HitSV/releases/tag/2.0.4\nSource-commit: %s\nSource-archive: %s\nLicense: GPL-3.0\nSource-scope: Release repository snapshot; core SVcalling source absent at pinned commit\nError-profile: 9mer_profile extracted from Source-archive demo/9mer_profile (official Hybrid -k profile)\n' "$URL" "$COMMIT" "$SOURCE_URL" > "$DEST/PROVENANCE.txt"
# Preserve the release repository snapshot and its GPLv3 license.
tar -xOf "$SOURCE_ARCHIVE" "HitSV-$COMMIT/LICENSE" > "$DEST/LICENSE"
