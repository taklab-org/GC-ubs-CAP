function [] = plot_manifold(a,rho2,nu)

%  Plots the graph of u(r) = sum_{k=0}^N a_k r^k
%
%  input: a = (a1,a2), where ai = (ai_0,ai_1,ai_2,...,ai_N)  (the Taylor coefficients)

figure

N=(length(a)-2)/2;

a1 = a(1:N+1); a2 = a(N+2:2*N+2);

t=(-nu:.001:nu);

m = length(t);

u1 = zeros(1,m); u2 = zeros(1,m);

for k = 1:m
    u1(k) = sum(a1.*t(k).^(0:N)');
    u2(k) = sum(a2.*t(k).^(0:N)');
end
    
plot(u1,u2,'Linewidth',3,'color',[1 0 0])
hold on
plot(rho2,0,'*','Linewidth',4,'color',[0 0 0])

set(gca,'FontSize',20)

axis tight
xlabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 30)

end

    