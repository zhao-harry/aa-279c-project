function K = getKnm(n,m)
    if n == 1
        K = 0;
    else
        K = ((n-1)^2 - m^2)/((2*n-1)*(2*n-3));
    end
end