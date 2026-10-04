#!/usr/bin/env bash
# Interleaved A/B of two hmmsearch builds: byte-identical outputs + speed. See ab.py --help.
exec python3 "$(dirname "${BASH_SOURCE[0]}")/ab.py" "$@"
