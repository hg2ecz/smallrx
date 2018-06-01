#!/bin/bash

dd if=/dev/urandom bs=1M count=30 of=testfile
pushd .; cd Rust; ./build.sh; popd
echo

for opt in '-Ofast' '-O3'; do
    # GCC
    gcc_compile="gcc   $opt -march=native rx.c -lm -o rx-gcc$opt";
    echo -e "\n$gcc_compile";   `$gcc_compile`;   cat testfile | time -p ./rx-gcc$opt   14080000 14070000 f > /dev/null;

    # CLANG
    clang_compile="clang $opt -march=native rx.c -lm -o rx-clang$opt";
    echo -e "\n$clang_compile"; `$clang_compile`; cat testfile | time -p ./rx-clang$opt 14080000 14070000 f > /dev/null;
done

# Rust
echo -e "\n./Rust/target/release/smallrx"; cat testfile | time -p ./Rust/target/release/smallrx 14080000 14070000 f > /dev/null;
