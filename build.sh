#!/usr/bin/env bash

rm -rf build
mkdir build
pushd build

export CC=clang
export CXX=clang++

cmake .. \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      -Wno-dev

make -j4

cp compile_commands.json ../

popd
