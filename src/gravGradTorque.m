function M = gravGradTorque(Ix,Iy,Iz,n,c)
    cx = c(1);
    cy = c(2);
    cz = c(3);
    M(1) = 3 * n^2 * (Iz - Iy) * cy * cz;
    M(2) = 3 * n^2 * (Ix - Iz) * cz * cx;
    M(3) = 3 * n^2 * (Iy - Ix) * cx * cy;
end