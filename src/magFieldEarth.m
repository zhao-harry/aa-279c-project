function [B_R,B_theta,B_phi] = magFieldEarth(R,phi,theta,RE)
    % NOTE: slides calls phi lambda (not sure if it's the same thing or not)

    % Make sure that R is normalized
    R = norm(R);

    % For g & h matrix (row = n, col = m + 1)
    g = [-30186 -2036 0 0 0; ...
        -1898 2997 1551 0 0; ...
        1299 -2144 1296 805 0; ...
        951 807 462 -393 235] * 1e-9; % T
    h = [0 5735 0 0 0; ...
        0 -2124 -37 0 0; ...
        0 -361 249 -253 0; ...
        0 148 -264 37 -307] * 1e-9; % T

    B_R = 0;
    B_theta = 0;
    B_phi = 0;
    for n = 1:4
        BR_temp = 0;
        BTheta_temp = 0;
        BPhi_temp = 0;

        for mInd = 1:n
            m = mInd - 1;

            P_nm = getPnm(theta,n,m);
            dPnm_dtheta = getdPdTheta(theta,n,m);

            BR_temp = BR_temp + ...
                (g(n,mInd) * cos(m * phi) + h(n,mInd) * sin(m * phi)) * ...
                P_nm;
            BTheta_temp = BTheta_temp + ...
                (g(n,mInd) * cos(m * phi)+ h(n,mInd) * sin(m * phi)) * ...
                dPnm_dtheta;
            BPhi_temp = BPhi_temp + ...
                (-g(n,mInd) * sin(m * phi) + h(n,mInd) * cos(m * phi)) * ...
                m * P_nm;
        end

        B_R = B_R + (RE / R)^(n + 2) * (n + 1) * BR_temp;
        B_theta = B_theta + (RE / R)^(n + 2) * BTheta_temp;
        B_phi = B_phi + (RE / R)^(n + 2) * BPhi_temp;
    end

    B_theta = -B_theta;
    B_phi = -1 / sin(theta) * B_phi;

end