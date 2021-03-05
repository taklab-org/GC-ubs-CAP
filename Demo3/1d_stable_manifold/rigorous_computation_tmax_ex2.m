function t_max = rigorous_computation_tmax_ex2(a,lambda,r0,nu)

N = (length(a)-2)/2; % a = (a1,a2) with ai = (ai_n)_{n=0}^N \in R^{N+1}

a1 = a(1:N+1); a2 = a(N+2:2*N+2);

A1 = cauchy_ext([a1 a1 a1 a1],4*N+1);
A2 = cauchy_ext([a2 a2],2*N+1);
A3 = cauchy_ext([a2 a2 a2 a2],4*N+1);
A4 = cauchy_ext([a1 a1 a1 a1 a2 a2],6*N+1);
A5 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1],8*N+1);

SIGMA = (-nu:.01:nu);
m = length(SIGMA);
t_max = intval(zeros(1,m));

for k=1:m-1
    sigma = SIGMA(k);
    
    na1 = norm(a1,1); na2 = norm(a2,1);
    
    n1 = intval(1:4*N)';
    sum1 = (sum((A1(2:end)./n1).*sigma.^n1)+4*(na1+r0)^3*r0)/lambda;
    
    n2 = intval(1:2*N)';
    sum2 = (sum((A2(2:end)./n2).*sigma.^n2)+2*(na2+r0)*r0)/lambda;
    
    n3 = intval(1:4*N)';
    sum3 = (sum((A3(2:end)./n3).*sigma.^n3)+4*(na2+r0)^3*r0)/lambda;
    
    n4 = intval(1:6*N)';
    error4 = 4*(na1+r0)^3*r0*(na2+r0)^2 + 2*(na1+r0)^4*(na2+r0)*r0;
    sum4 = (sum((A4(2:end)./n4).*sigma.^n4)+error4)/lambda;
    
    n5 = intval(1:8*N)';
    sum5 = (sum((A5(2:end)./n5).*sigma.^n5)+8*(na1+r0)^7*r0)/lambda;
    
    t_max(k) = sum1 + sum2 + sum3 + sum4 + sum5;
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

end