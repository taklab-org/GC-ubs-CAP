clf
clear
% figure
LW = 'linewidth'; lw = 1.6;
% n = 30;
% for r = 0.1:0.1:1
%   for th = pi/n:pi/n:2*pi
% p0 = [r*cos(th),r*sin(th)];
pb = [0.7328506362011802,0.5370700549804747];% pb
p_infs = [0.8861081289780320,0.6192579489210105];% p_{inf,s}
ps = [0.832267964033316,0.570370294663353]; % end point of stable manifold attached to p_{inf,s}
% 


opts = odeset('abstol',1e-18,'reltol',1e-18);

[~,y] = ode45(@vectorfield, [0,20], [0.88610810,0.61927], opts);
plot(y(:,1),y(:,2),'r',LW,lw), 
hold on 

% [~,y] = ode45(@vectorfield, [0,100], [-0.886109,0.6191], opts);
% plot(y(:,1),y(:,2),'r',LW,lw)

[~,y] = ode45(@vec_backward, [0,100], [0.886109,0.6191], opts);
plot(y(:,1),y(:,2),'r',LW,lw)

[~,y] = ode45(@vectorfield, [0,20], [0.88610813,0.61924], opts);
plot(y(:,1),y(:,2),'b',LW,lw)
% 
% [~,y] = ode45(@vec_backward, [0,20], [-0.886109,0.6191], opts);
% plot(y(:,1),y(:,2),'b',LW,lw)


for y0 = 0:0.05:1-eps
x0 = [0,y0];
% [~,y] = ode45(@vectorfield, [0,150], x0,opts);
% % 
% if y0 < 0
% plot(y(:,1),y(:,2),'b',LW,lw)
% else
% plot(y(:,1),y(:,2),'r',LW,lw)
% end
[~,y] = ode45(@vec_backward, [0,150], x0,opts);
if y0 < 0
plot(y(:,1),y(:,2),'b',LW,lw)
else
plot(y(:,1),y(:,2),'r',LW,lw)
end

% count = count + 1;
% end
%   end
end

p_sink = p_infs; p_infs(1) = -p_infs(1);
p_source = p_infs;
v_si = p_sink/norm(p_sink);
v_so = p_source/norm(p_source);

p_sink = [0.989136995894978, 0.206758557005181];
p_source = [-0.989136995894978, 0.206758557005181];

for j = 0.01:0.05:1
x0 = j*v_si;
% opts = odeset('abstol',1e-18,'reltol',1e-18);
[~,y] = ode45(@vectorfield, [0,150], x0,opts);
plot(y(:,1),y(:,2),'b',LW,lw)
% 
[~,y] = ode45(@vec_backward, [0,150], x0,opts);
plot(y(:,1),y(:,2),'b',LW,lw)
end

% for j = 0.01:0.05:1
% x0 = j*v_so;
% % opts = odeset('abstol',1e-18,'reltol',1e-18);
% [~,y] = ode45(@vectorfield, [0,150], x0,opts);
% plot(y(:,1),y(:,2),'b',LW,lw)
% % 
% [~,y] = ode45(@vec_backward, [0,150], x0,opts);
% plot(y(:,1),y(:,2),'b',LW,lw)
% end

p0 = [0, 0];% p0
% plot(p0(1),p0(2),'o','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
pb = [0.7328506362011802, 0.5370700549804747];% pb
% plot(pb(1),pb(2),'^','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(-pb(1),pb(2),'s','MarkerSize',20,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
p_infs = [0.8861081289780320, 0.6192579489210105];% p_{inf,s}
plot(p_infs(1),p_infs(2),'o','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(-p_infs(1),p_infs(2),'o','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(p_sink(1),p_sink(2),'s','MarkerSize',20,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
% plot(p_source(1),p_source(2),'^','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')
plot(pb(1),pb(2),'^','MarkerSize',16,'MarkerEdgeColor','k', 'MarkerFaceColor','k')

v = (p_infs-ps);
v_normal = [v(2),-v(1)]; 
v_normal = v_normal/norm(v_normal);
count = 0;
% p0 = [0.1,0.1];
resolution = 5e-3;
for t = 0:resolution:1
  % for t = -5e-3%:-5e-3:-1e-2
  p0 = ps + t*v_normal;
  if px4(p0) > 1
%     disp('hi')
    break
  end
  plot(p0(1),p0(2),'ko','MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k')
%   [~,y] = ode45(@vectorfield, [0,150],p0,opts);
%   plot(y(:,1),y(:,2),'g',LW,lw)
end



hold off

axis([0.68 1 0.43 0.67])
% axis equal

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