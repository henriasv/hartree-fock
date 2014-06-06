close all; clear all; clc

A = importdata('/scratch/hfdata/H2O_energies_431g.dat', ' ');

energies = A(:, 1);
rs = A(:, 2);
thetas = A(:, 3)*180/pi;

[val, ind] = min(energies)

F = TriScatteredInterp(rs,thetas,energies);

%for i=1:length(rs)-1
%    if (rs(i) == rs(i+1))
%        rs(i) = NaN;
%    end
%end
%for i=1:length(thetas)-1
%    if (thetas(i) == thetas(i+1))
%        thetas(i) = NaN;
%    end
%end
%rs = rs(rs==rs);
%thetas = thetas(thetas==thetas);

[X, Y] = meshgrid(sort(rs), sort(thetas));

Z = F(X, Y);

cm = [jet(128)];
pcolor(X, Y, Z);
shading interp

hold on
%caxis([-76.1 -75]);
colormap(cm);
colorbar

title('Energy as a function of nuclei configuration for H_2O with 3-21G basis', 'fontsize', 14)
set(gca, 'fontsize', 12);
xlabel('r [a.u.]');
ylabel('\theta [^\circ]');
h1 = scatter(rs(ind), thetas(ind), 'or', 'fill', 'displayname', 'Minimum value');
legend(h1,get(h1, 'displayname'))

%scatter3(rs, thetas, energies);