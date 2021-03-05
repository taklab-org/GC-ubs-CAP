function t_max = rigorous_computation_tmax_ex3_2d(a,scaling,fp,theta,r0)

N = round((-3+sqrt(25+8/3*length(a)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

i0 = (1:(N+1)*(N+2)/2-3); i1 = i0+(N+1)*(N+2)/2-3; i2 = i1+(N+1)*(N+2)/2-3;
a0 = a(i0); a1 = a(i1); a2 = a(i2);

load initial_data_ex3

if fp == 1   
    lambda = [tlambda1_1;tlambda1_2]; 
    v1 = scaling(1)*tv1_1; 
    v2 = scaling(2)*tv1_2;
    tx = ip1;
end

if fp == 2
    lambda = [tlambda2_1 tlambda2_2]; 
    v1 = scaling(1)*tv2_1; 
    v2 = scaling(2)*tv2_2;
    tx = ip2;
end
    
a0 = reshape_data_2d_manifold_vec2matrix(a0,tx(1),v1(1),v2(1));
a1 = reshape_data_2d_manifold_vec2matrix(a1,tx(2),v1(2),v2(2));
a2 = reshape_data_2d_manifold_vec2matrix(a2,tx(3),v1(3),v2(3));

a0_2 = cauchy2d_ext(a0,a0);
a1_2 = cauchy2d_ext(a1,a1);
a2_2 = cauchy2d_ext(a2,a2);

sum1 = a0_2 + a1_2 + a2_2;

sigma_gap = intval(min(abs(sup(lambda))));

theta1 = theta(1); theta2 = theta(2);
SUM = 0; 

for alpha1 = 0:2*N
    for alpha2 = 0:2*N
        if alpha1+alpha2>0
            SUM = SUM + sum1(alpha1+1,alpha2+1)*(theta1.^alpha1)*(theta2.^alpha2)/(alpha1*lambda(1)+alpha2*lambda(2));
        end
    end
end

tilde_delta = sup((2/sigma_gap)*(norm1(a0)+norm1(a1)+norm1(a2))*r0 + 3*r0^2/sigma_gap);

t_max = -SUM + infsup(-tilde_delta,tilde_delta);

end