#!/bin/bash
py=python3
nn_fit=../../../src/nn_fit.exe
$py write_gauss.py
$nn_fit nn_fit_1D.txt > nn_fit_1D.out

