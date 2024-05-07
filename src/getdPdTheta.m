function dPdTheta = getdPdTheta(theta, n, m)
    if n < m
        error("n >= m for Legendre functions")
    end
    if n == m && m == 0
        dPdTheta = 0;
    elseif n == m
        dPdTheta = sin(theta)*getdPdTheta(theta, n-1, n-1) + cos(theta)*getPnm(theta, n-1, n-1);
    else
        K = getKnm(n,m);
        if K == 0
            dPdTheta = cos(theta)*getdPdTheta(theta, n-1, m) - sin(theta)*getPnm(theta, n-1, m);
        else
            dPdTheta = cos(theta)*getdPdTheta(theta, n-1, m) - sin(theta)*getPnm(theta, n-1, m) - K*getdPdTheta(theta, n-2, m);
        end
    end
end