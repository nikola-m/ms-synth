#!/usr/bin/env python3
"""Thin wrapper kept for backwards compatibility.

Equivalent to the installed console script ``ms-synth-generate``.
Run ``ms-synth-generate --help`` for options.
"""
from ms_synth.cli import main

if __name__ == "__main__":
    main()
