function gradV = magFieldPotentialGrad(R, lambda, theta, RE, dR, dLambda, dTheta)
    % Find gradient of magnetic field potential

    % Find gradV in R
    fwd_step = magFieldPotential(R+dR/2, lambda, theta, RE);
    bwd_step = magFieldPotential(R-dR/2, lambda, theta, RE);
    gradV_R = (fwd_step - bwd_step)/(2*dR);

    % Find gradV in gamma
    fwd_step = magFieldPotential(R, lambda+dLambda/2, theta, RE);
    bwd_step = magFieldPotential(R, lambda-dLambda/2, theta, RE);
    gradV_Gamma = (fwd_step - bwd_step)/(2*dLambda);

    % Find gradV in lambda
    fwd_step = magFieldPotential(R, lambda, theta+dTheta/2, RE);
    bwd_step = magFieldPotential(R, lambda, theta-dTheta/2, RE);
    gradV_Lambda = (fwd_step - bwd_step)/(2*dTheta);

    gradV = [gradV_R, gradV_Gamma, gradV_Lambda];
end