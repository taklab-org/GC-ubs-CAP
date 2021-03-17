function t_max = rigorous_computation_tmax(a,rho1,rho2,r0,nu)

N = (length(a)-2)/2; % a = (a1,a2) with ai = (ai_n)_{n=0}^N \in R^{N+1}

a1 = a(1:N+1); a2 = a(N+2:2*N+2);

a1a1 = cauchy([a1 a1]);
q = [a2(2:end);0];

%%%% We now compute the Taylor expansion of R(u) = Q(u)/P1(u)^2

%%% We first refine the data for r using Newton

inv_a1a1 = invert(mid(a1a1));
r = mid(cauchy([q inv_a1a1]));

f = mid(f_solve_for_r(r,a1,q));
nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0; tol = 5e-16;
Df = mid(Df_solve_for_r(a1));

while (k<=100) && (nf > tol)
    r = r - Df\f;
    f = mid(f_solve_for_r(r,a1,q));
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

%%%%%%%%%%
%%% Y0 %%%
%%%%%%%%%%

inu = intval(nu);
iDf = Df_solve_for_r(a1);
A = inv(mid(iDf));

f_ext = f_solve_for_r_ext(r,a1,q);

Y0_1 = norma([A*f_ext(1:N+1);f_ext(N+2:end)/mid(a1a1(1))],inu);

norm_A = intval(max([sup(mnorma(A,inu)),sup(1/a1a1(1))]));

Y0_2 = norm_A*(2*norma(a1,inu)*norma(r,inu)+norma(r,inu)*r0+1/inu)*r0;

Y0 = sup(Y0_1+Y0_2);

disp(['Y0 = ',num2str(Y0)])

%%%%%%%%%%
%%% Z0 %%%
%%%%%%%%%%

iA = intval(A);
B = intval(eye(N+1) - iA*iDf);
Z0 = sup(mnorma(B,nu));

disp(['Z0 = ',num2str(Z0)])

%%%%%%%%%%
%%% Z1 %%%
%%%%%%%%%%

Z1_1 = norm_A*(2*r0*norma(a1,inu)+r0^2);
a1a1_ext = cauchy_ext([a1 a1]);
Z1_2 = (1/a1a1(1))*(norma([0;a1a1_ext(2:end)],inu)+2*norma(a1,inu)*r0+r0^2);

Z1 = sup(Z1_1+Z1_2);

disp(['Z1 = ',num2str(Z1)])

if sup(intval(Z0)+intval(Z1))<1
    disp('success with controling the r_n')
else
    disp('failure with controling the r_n')
    return
end

r_min = sup(intval(Y0)/(1-intval(Z0)-intval(Z1)));

lambda = -(rho2/2)*(rho2-rho1); % Stable eigenvalue

n = intval(0:N)';
SIGMA = (-nu:.01:nu);
m = length(SIGMA);
t_max = intval(zeros(1,m));
ir = intval(r);

for k=1:m-1
    %sigma = infsup(SIGMA(k),SIGMA(k+1));
    sigma = SIGMA(k);
    t_max(k) = (-1/lambda)*sum(ir./(n+1).*sigma.^(n+1))+r_min/abs(lambda)*infsup(-1,1);
end

plot(SIGMA(1:m-1),mid(t_max(1:m-1)),'Linewidth',3,'color',[0 0 0.8])

% %%%%% take data for plot3 of blow-up time
% u1 = zeros(1,m); u2 = zeros(1,m);
% for k = 1:m
%     u1(k) = sum(a1.*SIGMA(k).^(0:N)');
%     u2(k) = sum(a2.*SIGMA(k).^(0:N)');
% end
% x = u1(1:102); y = u2(1:102); z = mid(t_max(1:102));
% save('stable_manifold_xyz.mat','x','y','z')
% %%%%%

hold on
plot(0,0,'*','Linewidth',4,'color',[0 0 0])

set(gca,'FontSize',20)

axis tight
xlabel('$$\theta$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$t_{\max}(\theta)$$', 'Interpreter', 'latex', 'FontSize', 30)
grid on

end

function f = f_solve_for_r(r,a1,q)

a1a1r = cauchy([a1 a1 r]);
f = a1a1r - q;

end

function f = f_solve_for_r_ext(r,a1,q)

N = length(a1)-1 ;
a1a1r = cauchy_ext([a1 a1 r]);
f = a1a1r - [q;zeros(2*N,1)];

end

function Df = Df_solve_for_r(a1)

N = length(a1)-1 ;
Df = intval(zeros(N+1));

a1a1 = cauchy([a1 a1]);

for n = 0:N
    j = (0:n);
    k = n-j+1;
    Df(n+1,j+1) = a1a1(k);
end

end