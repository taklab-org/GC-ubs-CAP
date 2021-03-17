function plot_blowup_times(tmax)
fileID = fopen('data.bin');
A = fread(fileID,'double');
fclose(fileID);
xy = reshape(A,2,length(A)/2);

plot(xy(1,:),xy(2,:),'o','MarkerSize',6,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor',[.6 1 .6])
hold on
% % xlabel('$\gamma$','interpreter','latex','FontSize', 30)
% xlabel('Parametervalues','FontSize', 20)
% % ylabel('$$','interpreter', 'latex','FontSize', 30)
% ylabel('Blow-up time','FontSize', 20)
set(gca,'FontSize',15)
xlabel('Distance from the stable manifold','FontSize', 20)
% ylabel('$$','interpreter', 'latex','FontSize', 30)
ylabel('Blow-up time','FontSize', 20)

% load stable_manifold_xyz.mat
% 
plot(0,mid(tmax),'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
hold off
% 