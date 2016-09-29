#!/bin/bash

cd src
python3 setup.py build_ext --inplace
cd ..
ln -sf ./src/tdpost.*.so tdpost.so
