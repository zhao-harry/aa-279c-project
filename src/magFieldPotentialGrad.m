function gradV = magFieldPotentialGrad(R, lambda, theta, RE, dR, dLambda, dTheta)
    % Find gradient of magnetic field potential

    % Convert R to R_norm if needed
    if length(R) == 3
        R = norm(R);
    end

    % Find gradV in R
    fwd_step = magFieldPotential(R+dR, lambda, theta, RE);
    bwd_step = magFieldPotential(R-dR, lambda, theta, RE);
    gradV_R = (fwd_step - bwd_step)/(2*dR);

    % Find gradV in gamma
    fwd_step = magFieldPotential(R, lambda+dLambda, theta, RE);
    bwd_step = magFieldPotential(R, lambda-dLambda, theta, RE);
    gradV_Gamma = (fwd_step - bwd_step)/(2*dLambda);

    % Find gradV in lambda
    fwd_step = magFieldPotential(R, lambda, theta+dTheta, RE);
    bwd_step = magFieldPotential(R, lambda, theta-dTheta, RE);
    gradV_Lambda = (fwd_step - bwd_step)/(2*dTheta);

    gradV = [gradV_R, gradV_Gamma, gradV_Lambda];
end