function MJD = UT12MJD(UT1)
    % Input UT1 in YMD
    Y = UT1(1); M = UT1(2); D = UT1(3);
    if M <= 2
        y = Y - 1;
        m = M + 12;
    else
        y = Y;
        m = M;
    end
    B = floor(y / 4000) - floor(y / 100) + floor(y / 4);
    MJD = 365 * y - 679004 + floor(B) + floor(30.6001 * (m + 1)) + D;
end