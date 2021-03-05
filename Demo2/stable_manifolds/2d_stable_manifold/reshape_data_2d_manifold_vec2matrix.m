function ai_matrix = reshape_data_2d_manifold_vec2matrix(ai,txi,v1i,v2i)

N = round((-3+sqrt(25+8*length(ai)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

if isintval(ai(1,1)) == 1
    ai_matrix = intval(zeros(N+1));
else
    ai_matrix = zeros(N+1);
end
    
ai_matrix(1,1) = txi;
ai_matrix(2,1) = v1i;
ai_matrix(1,2) = v2i;

% We define the first two columns (including first order data)
ai_matrix(3:N+1,1) = ai(1:N-1); ai(1:N-1) = [];
ai_matrix(2:N,2) = ai(1:N-1); ai(1:N-1) = [];

% We define the rest of the upper triangular columns
for j = 2:N
    i = (0:N-j);
    ai_matrix(i+1,j+1) = ai(1:N-j+1);
    ai(1:N-j+1) = [];
end

end

