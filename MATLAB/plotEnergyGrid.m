close all; clear all; clc

rs = fread(fopen('/scratch/hfdata/H2O_rs_sto3g.bin', 'r'), 'double')';
thetas = fread(fopen('/scratch/hfdata/H2O_thetas_sto3g.bin', 'r'), 'double')';
N_t = length(thetas);
N_r = length(rs);
energies = fread(fopen('/scratch/hfdata/H2O_energies_sto3g.bin', 'r'), 'double');
energies = reshape(energies, N_t, N_r);


contour(energies, 100)
shading interp
colorbar