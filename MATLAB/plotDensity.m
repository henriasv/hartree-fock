filename = '/scratch/densityH2O_sym_1045deg431.bin';
a = fread(fopen(filename, 'r'), 'double');
size = length(a);
n = sqrt(size);
X = linspace(0, 6, n);
density = reshape(a, n, n);
figure(1);
Z = log10(density);
Z = Z(1:end, 1:end);
%pcolor(Z)
cm = [jet(64);gray(0)];
shading flat
hold all
contour(X, X, Z, 100, 'linewidth', 1)
axis square
colormap(cm);
title('Logarithm of particle density for H_2O with 4-31G basis', 'fontsize', 14)
set(gca, 'fontsize', 12);
xlabel('x [a.u.]');
ylabel('y [a.u.]');
colorbar;
