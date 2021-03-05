function [] = plot_manifold(a,tx,nu)

%  Plots the graph of u(r) = sum_{k=0}^N a_k r^k
%
%  input: a = (a0,a1,a2), where ai = (ai_0,ai_1,ai_2,...,ai_N)  (the Taylor coefficients)

%figure

N=(length(a)-3)/3;

a0 = a(1:N+1); a1 = a(N+2:2*N+2); a2 = a(2*N+3:3*N+3);

t=(-nu:.001:nu);

m = length(t);

u0 = zeros(1,m); u1 = zeros(1,m); u2 = zeros(1,m);

for k = 1:m
    u0(k) = sum(a0.*t(k).^(0:N)');
    u1(k) = sum(a1.*t(k).^(0:N)');
    u2(k) = sum(a2.*t(k).^(0:N)');
end
    
plot3(u0,u1,u2,'Linewidth',3,'color',[1 0 0])
hold on
plot3(tx(1),tx(2),tx(3),'*','Linewidth',4,'color',[0 0 0])

set(gca,'FontSize',20)

axis tight
xlabel('$$x_0$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 30)
zlabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 30)

end

    