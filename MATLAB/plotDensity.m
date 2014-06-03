filename = '/scratch/densityH2.bin';
a = fread(fopen(filename, 'r'), 'double');
size = length(a);
n = sqrt(size);
X = linspace(0, 4, n);
density = reshape(a, n, n);
figure(1);
Z = log10(density);
Z = Z(1:end, 1:end);
%pcolor(Z)
cm = [jet(64);gray(0)];
shading flat
hold all
contour(X, X, Z, 100, 'linewidth', 2)
axis square
colormap(cm);
title('Logarithm of particle density for H_2 with STO-3G basis', 'fontsize', 14)
set(gca, 'fontsize', 12);
colorbar;