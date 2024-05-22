function rNoise = addVectorNoise(rNorm, angLim)
    theta = acos(rNorm(3));
    phi = sign(rNorm(2)) * acos(rNorm(1) / norm(rNorm(1:2)));

    theta = theta + angLim*randn;
    phi = phi + angLim*randn;

    rNoise = [sin(theta)*cos(phi);
                   sin(theta)*sin(phi);
                   cos(theta)];
end