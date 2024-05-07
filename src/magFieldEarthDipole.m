function B = magFieldEarthDipole(R,mE,rE3_B0)
    R_norm = norm(R);
    B = -rE3_B0 / R_norm^3 * (3 * dot(mE,R) * R - mE);
end