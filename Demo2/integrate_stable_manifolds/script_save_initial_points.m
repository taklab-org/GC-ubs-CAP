
% load p1.mat
% fileID = fopen('initial_points_p1.bin','w');

load p2.mat
fileID = fopen('initial_points_p2.bin','w');

pp_inf = inf(pp); pp_inf = pp_inf(:);
pp_sup = sup(pp); pp_sup = pp_sup(:);

N = length(pp_inf);

init_data = zeros(2*N,1);

for i = 1:N
  init_data(2*(i-1)+1)   = pp_inf(i);
  init_data(2*i) = pp_sup(i);
end
% 
fwrite(fileID,init_data,'double');
fclose(fileID);