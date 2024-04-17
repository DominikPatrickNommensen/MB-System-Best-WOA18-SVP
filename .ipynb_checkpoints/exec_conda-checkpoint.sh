#!/bin/bash

# Extract arguments
args="$@"

# Run with Conda
conda run -n svptool --live-stream python ../../MB-System-Best-WOA18-SVP/mbbestsvp23.py $args
