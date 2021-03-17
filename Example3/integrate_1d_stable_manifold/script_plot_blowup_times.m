% function plot_blowup_times(tmax)
clf
fileID = fopen('data.bin');
A = fread(fileID,'double');
fclose(fileID);
xy = reshape(A,2,length(A)/2);

plot(xy(1,:),xy(2,:),'o','MarkerSize',6,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
hold on
% xlabel('$\gamma$','interpreter','latex','FontSize', 30)
set(gca,'FontSize',15)
xlabel('Distance from the stable manifold','FontSize', 20)
% ylabel('$$','interpreter', 'latex','FontSize', 30)
ylabel('Blow-up time','FontSize', 20)

% load stable_manifold_xyz.mat
% % % 
% plot(0,mid(tmax),'o','MarkerSize',10,...
%     'MarkerEdgeColor','red',...
%     'MarkerFaceColor',[1 .6 .6])
% hold off
% % 3.0845