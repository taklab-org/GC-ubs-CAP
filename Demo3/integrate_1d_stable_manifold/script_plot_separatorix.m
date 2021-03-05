close all
% clf
% clear
% figure
open stable_manifold.fig
hold on

LW = 'linewidth'; lw = 2;
% n = 30;
% for r = 0.1:0.1:1
%   for th = pi/n:pi/n:2*pi
% p0 = [r*cos(th),r*sin(th)];
% x0 = [0.7328506362011802,0.5370700549804747];% pb
p_infs = [0.8861081289780320,0.6192579489210105];% p_{inf,s}
ps = [0.832267964033316,0.570370294663353]; % end point of stable manifold attached to p_{inf,s}
% 

opts = odeset('abstol',1e-18,'reltol',1e-18);
[~,y] = ode45(@vec_backward, [0,150],ps,opts);
plot(y(:,1),y(:,2),':k','Linewidth',3)

v = (p_infs-ps);
v_normal = [v(2),-v(1)]; 
v_normal = v_normal/norm(v_normal);
% count = 0;
% p0 = [0.1,0.1];
% resolution = 5e-3;
% for t = 0:resolution:1
for t = [5e-3,-5e-3]%:-5e-3:-1e-2
p0 = ps + t*v_normal;
if px4(p0) > 1
  break
end
% x = p0;
% p0(1) = p0(1)/(1-px4(x));
% p0(2) = p0(2)/(1-px4(x))^2;

% eps = 1e-16;
opts = odeset('abstol',1e-18,'reltol',1e-18);
[~,y] = ode45(@vectorfield, [0,150], p0,opts);
% [t,y] = ode45(@vec_original, [0,10], p0, opts);
if t > 0
plot(y(:,1),y(:,2),'b',LW,lw)
else
plot(y(:,1),y(:,2),'r',LW,lw)
end
hold on
% count = count + 1;
end


ps = [0.320682556228925, 0.319896988077054];

opts = odeset('abstol',1e-18,'reltol',1e-18);
[~,y] = ode45(@vec_backward, [0,150],ps,opts);
plot(y(:,1),y(:,2),':k','Linewidth',3)

v = -ps;
v_normal = [v(2),-v(1)]; 
v_normal = v_normal/norm(v_normal);
% count = 0;
% p0 = [0.1,0.1];
resolution = 5e-3;
% for t = 0:resolution:1
for t = [5e-3,-5e-3]%:-5e-3:-1e-2
p0 = ps + t*v_normal;
if px4(p0) > 1
  break
end
opts = odeset('abstol',1e-18,'reltol',1e-18);
[~,y] = ode45(@vectorfield, [0,150], p0,opts);
% [t,y] = ode45(@vec_original, [0,10], p0, opts);
if t > 0
plot(y(:,1),y(:,2),'r', LW, lw)
else
plot(y(:,1),y(:,2),'b',LW,lw)
end
% hold on
% count = count + 1;
end
% 
% p_sink = p_infs; p_infs(1) = -p_infs(1);
% p_source = p_infs;
% v_si = p_sink/norm(p_sink);
% v_so = p_source/norm(p_source);
% 
% p_sink = [0.989136995894978, 0.206758557005181];
% p_source = [-0.989136995894978, 0.206758557005181];

% for j = 0.01:0.05:1
% x0 = j*v_si;
% % opts = odeset('abstol',1e-18,'reltol',1e-18);
% [t,y] = ode45(@vectorfield, [0,150], x0,opts);
% plot(y(:,1),y(:,2),'b',LW,lw)
% % 
% [t,y] = ode45(@vec_backward, [0,150], x0,opts);
% plot(y(:,1),y(:,2),'b',LW,lw)
% end

% for j = 0.01:0.05:1
% x0 = j*v_so;
% % opts = odeset('abstol',1e-18,'reltol',1e-18);
% [t,y] = ode45(@vectorfield, [0,150], x0,opts);
% plot(y(:,1),y(:,2),'b',LW,lw)
% % 
% [t,y] = ode45(@vec_backward, [0,150], x0,opts);
% plot(y(:,1),y(:,2),'b',LW,lw)
% end

p0 = [0, 0];% p0
plot(p0(1),p0(2),'o','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')

pb = [0.7328506362011802, 0.5370700549804747];% pb
plot(pb(1),pb(2),'^','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
plot(-pb(1),pb(2),'s','MarkerSize',20,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
p_infs = [0.8861081289780320, 0.6192579489210105];% p_{inf,s}
plot(p_infs(1),p_infs(2),'o','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(-p_infs(1),p_infs(2),'o','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(p_sink(1),p_sink(2),'s','MarkerSize',20,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
p_sink = [0.989136995894978, 0.206758557005181];
plot(p_sink(1),p_sink(2),'s','MarkerSize',20,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(p_source(1),p_source(2),'^','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% 
% hold off


% separatrix: Figure 14
% x=[0.63,0.58];y=[0.5,0.4];annotation('textarrow',x,y);
% x=[0.6,0.52];y=[0.55,0.47];annotation('textarrow',x,y);
% x=[0.85,0.88];y=[0.65,0.53];annotation('textarrow',x,y);
% x=[0.82,0.77];y=[0.73,0.83];annotation('textarrow',x,y);

% % separatrix zoom1: Figure 15(a)
% grid on
% axis([0.26,0.38,0.26,0.38])
% axis square
% xticks([0.26:0.02:0.38])

% separatrix zoom2: Figure 15(b)
% grid on
% axis([0.8,0.92,0.52,0.64])
% xticks([0.8:0.02:0.92])
% axis square

set(gca,'FontSize',20)

xlabel('$$x_1$$', 'Interpreter', 'latex', 'FontSize', 30)
ylabel('$$x_2$$', 'Interpreter', 'latex', 'FontSize', 30)

function y = px4(x)
y = x(1)^4 + x(2)^2;
end

function y = F(x)
y = 0.25*(1. + 3*px4(x));
end

function y = G(x)
x1 = x(1); x2 = x(2);
y = x1.^3.*(x1.^2-x2) + 0.5*x2.*(x1.^3/3 - x1.*(1-px4(x)).^2);
end

function dx = vectorfield(t,x)
x1 = x(1); x2 = x(2); %dx = zeros(2,1);
dx = [(x1.^2 - x2).*F(x) - x1.*G(x);...
  (x1.^3/3 - x1.*(1-px4(x)).^2).*F(x) - 2*x2.*G(x)];
end

function dx = vec_backward(t,x)
dx = -vectorfield(t,x);
end

function dx = vec_original(t,x)
dx = [x(1).^2 - x(2); x(1).^3/3 - x(1)];
end