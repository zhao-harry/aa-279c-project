function [alpha, AE_GT]= getAlpha(qEstimated, A_ECI2P,rECI,vECI,control)
    h = cross(rECI,vECI);
    radial = rECI / norm(rECI);
    normal = h / norm(h);
    tangential = cross(normal,radial);
    A_Target = [-radial -normal -tangential]';
    AE_GT = A_ECI2P * A_Target';

    AE = q2A(qEstimated) * A_Target';
    % AE = AE_GT;

    if control.useLinearModel == true
        alpha = [AE(2,3); AE(3,1); AE(1,2)];
    else
        alpha = [(AE(2,3) - AE(3,2)) / 2; ...
            (AE(3,1) - AE(1,3)) / 2; ...
            (AE(1,2) - AE(2,1)) / 2];
    end
end