#!/bin/sh

RUSTFLAGS="-C target-cpu=native -C -Ofast" cargo build --release
