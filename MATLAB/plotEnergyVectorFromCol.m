
A = importdata('/scratch/hfdata/angularH2O_431g.dat', ' ');
B = importdata('/scratch/hfdata/angularH2O_321g.dat', ' ');

energies = A(:, 1);
rs = A(:, 2);
thetas = A(:, 3)*180/pi;
[val, ind] = min(energies);
[val2, ind2] = min(B(:, 1));

figure()
plot(thetas, energies);
hold all 
plot(B(:,3)*180/pi, B(:,1))
title('Angular dependence of energy for H_2O', 'fontsize', 14)
set(gca, 'fontsize', 12);
legend('4-31G', '3-21G');
scatter(thetas(ind), val, 'or', 'fill');
scatter(B(ind2, 3)*180/pi, val2, 'or', 'fill');

xlabel('\theta');
ylabel('E');

%%
C = importdata('/scratch/hfdata/distanceH2O_431g.dat', ' ');
D = importdata('/scratch/hfdata/distanceH2O_321g.dat', ' ');

[val, ind] = min(C(:, 1));
[val2, ind2] = min(D(:, 1));

figure()
plot(C(:, 2), C(:, 1));
hold all
plot(D(:, 2), D(:, 1));
ylim([-76 -75])

title('Bond distance dependence of energy for H_2O', 'fontsize', 14)
set(gca, 'fontsize', 12)
xlabel('r')
ylabel('E')
legend('4-31G', '3-21G');
scatter(C(ind, 2), val, 'or', 'fill');
scatter(D(ind2, 2), val2, 'or', 'fill');

%%
E = importdata('/scratch/hfdata/H2O_energies_431g.dat');
