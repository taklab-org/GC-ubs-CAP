close all
clear

fileID = fopen('data.bin');
A = fread(fileID,'double');
fclose(fileID);
xyz = reshape(A,9,length(A)/9);
plot3(xyz(3,:),xyz(6,:),xyz(9,:),'linewidth',3), hold on

load stable_manifold_xyz.mat

plot3(x,y,z,'linewidth',3), hold on

plot3(2,0,0,'ro','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
text(2,0,.5,'$$p_2$$','interpreter','latex','FontSize', 30)
xlabel('$x_1$','interpreter','latex','FontSize', 30)
ylabel('$x_2$','interpreter', 'latex','FontSize', 30)
zlabel('Blow-up time','interpreter', 'latex','FontSize', 20)
title('Solution profile $(x_1,x_2)$ and blow-up time','interpreter', 'latex','FontSize', 20)
% grid on
% hold off
points = [1,600, 2200, 3600, 5123];
for i = 1:5
index = points(i);%= [1,600, 2200, 3600, 5123]
plot3(xyz(3,index),xyz(6,index),xyz(9,index),'bo','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
x = infsup(xyz(1,index),xyz(2,index));
% disp(['x = ', dispintval(x)])
y = infsup(xyz(4,index),xyz(5,index));
% disp(['y = ', dispintval(x)])
tmax = infsup(xyz(7,index),xyz(8,index));
% disp(['tmax = ', dispintval(x)])
str1 = ['$P_',num2str(i),'$ & '];
str2 = ['$',dispintval(x),'$ & '];
str3 = ['$',dispintval(y),'$ & '];
str4 = ['$',dispintval(tmax),'$\\'];
disp([str1,str2,str3,str4])
end
hold off
grid on
set(gcf,'Position',[0,0,600,350])
text(xyz(3,points(1)),xyz(6,points(1)),xyz(9,points(1))+0.3,'$$P_1$$','interpreter','latex','FontSize', 25)
text(xyz(3,points(2)),xyz(6,points(2)),xyz(9,points(2))+0.3,'$$P_2$$','interpreter','latex','FontSize', 25)
text(xyz(3,points(3))-0.018,xyz(6,points(3))+0.007,xyz(9,points(3))+0.3,'$$P_3$$','interpreter','latex','FontSize', 25)
text(xyz(3,points(4))-0.012,xyz(6,points(4))-0.05,xyz(9,points(4))+0.4,'$$P_4$$','interpreter','latex','FontSize', 25)
text(xyz(3,points(5))-0.015,xyz(6,points(5)),xyz(9,points(5))+0.45,'$$P_5$$','interpreter','latex','FontSize', 25)

figure
plot3(xyz(9,end)-xyz(9,:),xyz(3,:),1./xyz(6,:),'linewidth',3),hold on
view([65 20])
grid on

xlabel('$t$','interpreter','latex','FontSize', 30)
ylabel('$\beta(t)$','interpreter', 'latex','FontSize', 30)
zlabel('$v(t)$','interpreter', 'latex','FontSize', 30)


for i = 1:5
index = points(i);%= [1,600, 2200, 3600, 5123]
plot3(xyz(9,end)-xyz(9,index),xyz(3,index),1./xyz(6,index),'bo','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
% x = infsup(xyz(1,index),xyz(2,index));
% disp(['x = ', dispintval(x)])
% y = infsup(xyz(4,index),xyz(5,index));
% disp(['y = ', dispintval(x)])
% tmax = infsup(xyz(7,index),xyz(8,index));
% disp(['tmax = ', dispintval(x)])
% str1 = ['$P_',num2str(i),'$ & '];
% str2 = ['$',dispintval(x),'$ & '];
% str3 = ['$',dispintval(y),'$ & '];
% str4 = ['$',dispintval(tmax),'$\\'];
% disp([str1,str2,str3,str4])
end

% set(gcf,'Position',[0,0,600,350])
title('Blow-up profile','interpreter', 'latex','FontSize', 25)

text(xyz(9,end)-xyz(9,points(1))+0.3, xyz(3,points(1)), 1./xyz(6,points(1)),'$$Q_1$$','interpreter','latex','FontSize', 25)
text(xyz(9,end)-xyz(9,points(2))+0.3, xyz(3,points(2)), 1./xyz(6,points(2)),'$$Q_2$$','interpreter','latex','FontSize', 25)
text(xyz(9,end)-xyz(9,points(3)), xyz(3,points(3))+0.009, 1./xyz(6,points(3)),'$$Q_3$$','interpreter','latex','FontSize', 25)
text(xyz(9,end)-xyz(9,points(4))+0.4, xyz(3,points(4))-0.012, 1./xyz(6,points(4))+2.5,'$$Q_4$$','interpreter','latex','FontSize', 25)
text(xyz(9,end)-xyz(9,points(5))+0.45, xyz(3,points(5))-0.028, 1./xyz(6,points(5))+1,'$$Q_5$$','interpreter','latex','FontSize', 25)

