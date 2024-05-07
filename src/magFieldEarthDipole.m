function B = magFieldEarthDipole(R, RE, mE, B0)
    R_norm = norm(R);
    B = - RE^3*B0/R_norm^3 * (3*dot(mE,R)*R - mE);
end