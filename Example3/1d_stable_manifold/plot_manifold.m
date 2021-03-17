function [] = plot_manifold(a,tx,nu,plus)

%  Plots the graph of u(r) = sum_{k=0}^N a_k r^k
%
%  input: a = (a1,a2), where ai = (ai_0,ai_1,ai_2,...,ai_N)  (the Taylor coefficients)

%figure

N=(length(a)-2)/2;

a1 = a(1:N+1); a2 = a(N+2:2*N+2);


if plus == 1
  t=(0:.001:nu);
elseif plus == -1
  t=(-nu:.001:0);
else
  t=(-nu:.001:nu);
end

m = length(t);

u1 = zeros(1,m); u2 = zeros(1,m);

% [sum(a1.*(t(1).^(0:N)')),sum(a2.*(t(1).^(0:N)'))]
% [sum(a1.*(t(m).^(0:N)')),sum(a2.*(t(m).^(0:N)'))]

for k = 1:m
    u1(k) = sum(a1.*t(k).^(0:N)');
    u2(k) = sum(a2.*t(k).^(0:N)');
end
    
if plus==0
  plot(u1,u2,'Linewidth',3,'color',[0.8 0 0])
else
  plot(u1,u2,'k:','Linewidth',3)
end

hold on
plot(tx(1),tx(2),'*','Linewidth',5,'color',[0 0 0])

set(gca,'FontSize',20)

axis tight
xlabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 30)

end

    