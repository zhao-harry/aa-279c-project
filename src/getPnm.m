function P = getPnm(theta,n,m)
    if n == 0 && m == 0
        P = 1;
    elseif n == m
        P = sin(theta) * getPnm(theta,n-1,n-1);
    else
        K = getKnm(n,m);
        if K == 0
            P = cos(theta) * getPnm(theta,n-1,m);
        else
            P = cos(theta) * getPnm(theta,n-1,m) - K * getPnm(theta,n-2,m);
        end
    end
end