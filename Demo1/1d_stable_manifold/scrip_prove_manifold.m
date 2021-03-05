close all
clear
load sol

nu = 1;

[I,success] = rad_poly(a,par,nu);

disp(['||P(theta) - P^(N)(theta)|| < = r_min = ',num2str(I(1))])