% clf
% clear
% figure
LW = 'linewidth'; lw = 1.4;

open 2d_stable_manifold_ex2.fig

hold on

a = 0.3; c = 0.7; delta = 9; w = 0.02;
p0 = [0.9333789,0.3588924,0];
p1 = [0.7180928,0.6959473,0];
p2 = [0.9985628,-0.0535924,0];
p_source = [0.7071051816183367,0.001504037399468,-0.001504037399468];

plot_trajectory([0.876149668193990,0.305009506850414,-0.000441867190016])

load p1.mat
for i = 1:length(pp)
  plot_trajectory(mid(pp(:,i)))
end

load p2.mat
for i = 1:length(pp)
  plot_trajectory(mid(pp(:,i)))
end

% plot_trajectory(mid(p2_1))
% plot_trajectory(mid(p2_2))
% plot_trajectory(mid(p2_3))
%     end
%   end
% end

% plot3(pp(1,1).mid,pp(2,1).mid,pp(3,1).mid,'o','MarkerSize',10,...
%     'MarkerEdgeColor','green',...
%     'MarkerFaceColor',[.6 1 .6])
%   
  
% plot3(p0(1),p0(2),p0(3),'o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor',[.6 .6 1])
% plot3(p1(1),p1(2),p1(3),'o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor',[.6 .6 1])
% plot3(p2(1),p2(2),p2(3),'o','MarkerSize',10,...
%     'MarkerEdgeColor','blue',...
%     'MarkerFaceColor',[.6 .6 1])
% 
% plot3(p_source(1),p_source(2),p_source(3),'o','MarkerSize',10,...
%     'MarkerEdgeColor','green',...
%     'MarkerFaceColor',[.6 1 .6])
%   
% xlabel('$x$','interpreter','latex','FontSize', 30)
% ylabel('$y$','interpreter', 'latex','FontSize', 30)
% zlabel('$z$','interpreter', 'latex','FontSize', 30)
view([-100,60])


function dx = vectorfield(t,x)
dx = g_vector_field(x);
end

function dx = vec_backward(t,x)
dx = -vectorfield(t,x);
end

function dx = vec_original(t,x,a,c,delta,w)
x1 = x(1); x2 = x(2); x3 = x(3);
dx = [x1.*(x1.^2 - 1);...
  (x1.^2).*(x2 + x3);...
  (x1.^2).*x3 + (c*(x1.^2).*x3 - x2.*(x2 - a*x1).*(x1 - x2) + w*x1.^3)/delta];
end

function plot_trajectory(x0)
opts = odeset('abstol',1e-18,'reltol',1e-18);
% [t,y] = ode45(@vectorfield, [0,200], x0,opts);
% [t,y] = ode45(@vectorfield, [0,36.588099706], p0,opts);
[~,y] = ode45(@vec_backward, [0,100], x0, opts);
% p0(1) = p0(1)/(1-px4(p0));
% p0(2) = p0(2)/(1-px4(p0))^2;
% [t,y] = ode45(@vec_original, [0,7.1], p0,opts);
% plot3(x0(1),x0(2),x0(3),'o','MarkerSize',6,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6]), hold on
plot3(y(:,1),y(:,2),y(:,3),'k','linewidth',1.6), hold on
plot3(x0(1),x0(2),x0(3),'o','MarkerSize',10,...
'MarkerEdgeColor','red',...
'MarkerFaceColor',[1 .6 .6])
end
