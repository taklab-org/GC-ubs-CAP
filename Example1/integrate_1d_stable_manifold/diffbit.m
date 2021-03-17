function n = diffbit(a,b)

n = 1;

while strncmpi(a,b,n)
    n = n+1;
end
% n = n-1;
end