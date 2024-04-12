function yECI = oe2eci(a,e,i,O,w,nu)
    i = deg2rad(i);
    O = deg2rad(O);
    w = deg2rad(w);
    nu = deg2rad(nu);
    p = a * (1 - e^2);
    r = p / (1 + e * cos(nu));
    rPQW = [r * cos(nu); r * sin(nu); 0];
    vPQW = sqrt(3.986 * 10^5 / p) * [-sin(nu); e + cos(nu); 0];
    Rzw = [cos(-w), sin(-w), 0;...
        -sin(-w), cos(-w), 0;...
        0, 0, 1];
    Rxi = [1, 0, 0;...
        0, cos(-i), sin(-i);...
        0, -sin(-i), cos(-i)];
    RzO = [cos(-O), sin(-O), 0;...
        -sin(-O), cos(-O), 0;...
        0, 0, 1];
    rECI = RzO * Rxi * Rzw * rPQW;
    vECI = RzO * Rxi * Rzw * vPQW;
    yECI = [rECI; vECI];
end