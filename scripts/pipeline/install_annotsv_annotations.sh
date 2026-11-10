#!/usr/bin/env bash
# Download the AnnotSV human annotations and Exomiser data into
# <DEST>/AnnotSV_annotations, for use with AnnotSV's -annotationsDir option.
# These data are not included in the AnnotSV Conda package or container.
#
# Usage: install_annotsv_annotations.sh [VERSION] [DEST]
#   VERSION  AnnotSV release, without the leading "v" (default: 3.5.3, the
#            version of the default AnnotSV container)
#   DEST     directory that will contain AnnotSV_annotations (default: .)
#
# The download is done by the install targets of AnnotSV's own Makefile, as
# documented in the AnnotSV README. Requires git, make, tclsh, curl, tar and
# unzip.
set -euo pipefail

version=${1:-3.5.3}
dest=$(cd "${2:-.}" && pwd)
target="$dest/AnnotSV_annotations"

if [[ -d "$target/Annotations_Human" ]]; then
    echo "AnnotSV annotations already installed in $target"
    exit 0
fi

for tool in git make tclsh curl tar unzip; do
    command -v "$tool" >/dev/null || { echo "Error: $tool is required" >&2; exit 1; }
done

work=$(mktemp -d "$dest/.annotsv_install.XXXXXX")
trap 'rm -rf "$work"' EXIT

echo "Fetching AnnotSV v$version ..."
git clone --quiet --depth 1 --branch "v$version" https://github.com/lgmgeo/AnnotSV.git "$work/src"

# Install inside the clone (PREFIX=.), the layout AnnotSV's Makefile supports for
# in-place use; only the annotation directories are kept. The "install" target
# is not needed: it pip-installs variantconvert, which HitFlow does not use.
(cd "$work/src" && make PREFIX=. install-human-annotation)

mkdir -p "$target"
for name in Annotations_Human Annotations_Exomiser; do
    mv "$work/src/share/AnnotSV/$name" "$target/"
done
echo "AnnotSV $version annotations installed in $target"
