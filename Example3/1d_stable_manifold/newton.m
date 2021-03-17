function a = newton(a,par)

tol = 5e-12; %% tolerance for Newton's method

f = f_stable_manifold(a,par);

nf = norm(f,1);
display(['At the beginning ||f|| = ',num2str(nf)])

k=0;

while (k<=100) && (nf > tol)
    %Df = mid(finite_diff_DF(a,par));
    Df = mid(Df_stable_manifold(a,par));
    a = a - Df\f;
    f = f_stable_manifold(a,par);
    nf = norm(f);
    display(['||f|| = ',num2str(nf),', ||Df^(-1)|| = ',num2str(norm(inv(Df),1))])
    k = k+1;                            
end

scaling = par(1); fp = par(2);

if fp == 0
    
    %%%%%%%%%%%%%%%
    %%%% (0,0) %%%%
    %%%%%%%%%%%%%%%
    
    tx = [0;0]; % Fixed point
    v = scaling*[1;1]; % Stable eigenvector

end

if fp == 1
   
    %%%%%%%%%%%%%%%%%%%%%
    %%%% p_{infty,s} %%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Fixed point (+ rigorous error bound for the proof)
    tx = [0.8861081289780320;0.6192579489210105];
    
    % Stable eigenvector (+ rigorous error bound for the proof)
    v = scaling*[-0.581870201791746;-0.813281666009282];
    
end

if fp == 2
   
    %%%%%%%%%%%%%
    %%%% p_b %%%%
    %%%%%%%%%%%%%
    
    % Fixed point (+ rigorous error bound for the proof)
    tx = [0.732850636201180;0.537070054980475];
    
    % Stable eigenvector (+ rigorous error bound for the proof)
    v = scaling*[-0.981046226514401;-0.193773840963773];
    
end

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}
a1 = a(1:N-1); a2 = a(N:2*N-2);
a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];
plot_manifold([a1;a2],tx,1);

end