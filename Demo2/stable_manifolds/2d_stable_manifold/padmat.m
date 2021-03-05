function [ai_ext] = padmat(ai,m)

N = length(ai) - 1;

if isintval(ai(1,1)) == 1
    ai_ext = intval(zeros(N+1+m));
else
    ai_ext = zeros(N+1+m);
end

ai_ext(1:N+1,1:N+1) = ai;

end

