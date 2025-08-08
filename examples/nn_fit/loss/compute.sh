#!/bin/bash
../../../src/nn_fit.exe nn_fit_1D_mse.txt
mv nn_fit_train.dat nn_mse_train.dat
mv nn_fit_val.dat nn_mse_val.dat
mv nn_fit_test.dat nn_mse_test.dat

../../../src/nn_fit.exe nn_fit_1D_huber.txt
mv nn_fit_train.dat nn_huber_train.dat
mv nn_fit_val.dat nn_huber_val.dat
mv nn_fit_test.dat nn_huber_test.dat

../../../src/nn_fit.exe nn_fit_1D_asinh.txt
mv nn_fit_train.dat nn_asinh_train.dat
mv nn_fit_val.dat nn_asinh_val.dat
mv nn_fit_test.dat nn_asinh_test.dat

