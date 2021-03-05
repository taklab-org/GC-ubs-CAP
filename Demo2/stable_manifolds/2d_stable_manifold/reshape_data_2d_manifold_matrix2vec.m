function ai = reshape_data_2d_manifold_matrix2vec(ai_matrix)

N = length(ai_matrix) - 1; % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

% if isintval(ai_matrix(1,1)) == 1
%     ai = intval(zeros((N+1)*(N+2)/2-3,1));
% else
%     ai = zeros((N+1)*(N+2)/2-3,1);
% end

ai(1:N-1,1) = ai_matrix(3:N+1,1);
ai(N:2*N-2,1) = ai_matrix(2:N,2);

for j = 2:N
    i = (0:N-j);
    ai = [ai;ai_matrix(i+1,j+1)];
end
 
% % We define the rest of the upper triangular columns
% for j = 2:N
%     i = (0:N-j);
%     ai_matrix(i+1,j+1) = ai(1:N-j+1);
%     ai(1:N-j+1) = [];
% end

end

