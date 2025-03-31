#!/usr/bin/env bash

echo "Removing spaces from: "
echo $1
one=$1
mv $one nolines_$one
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' nolines_$one > $one
rm nolines_$one
