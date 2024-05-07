function V = magFieldPotential(R, lambda, theta, RE)
    % Calculate the magnetic field potential
    % Inputs:
    % - R: position vector of satellite
    % - lambda: longitude of satellite
    % - theta: colatitude of satellite
    % - RE: radius of Earth
    
    % for g & h matrix (row = n, col = m + 1)
    g = [-30168, -2036, 0, 0, 0;
            -1898, 2997, 1551, 0, 0;
            1299, -2144, 1298, 805, 0;
            951, 807, 462, -393, 238] * 1e-9; % T
    h = [0, 5735, 0, 0, 0;
            0, -2124, -37, 0, 0;
            0, -361, 249, -253, 0;
            0, 148, -264, 37, -307] * 1e-9; % T

    % Find P matrix
    V_sum = 0;
    for n = 1:4
        innerTerm = 0;
        for mInd = 1:n
            m = mInd - 1;
            delta_0m = (m == 0); %Q: is 0 1 here because matlab indexing

            % calculate 2n-1
            fact_2nMin1 = 1;
            for k = 1:n
                fact_2nMin1 = fact_2nMin1 * (2*k-1);
            end

            P = (sqrt((2-delta_0m)*factorial(n-m)/factorial(n+m))*fact_2nMin1/factorial(n-m))*(sin(theta))^m ...
                    * ((cos(theta))^(n-m) - (n-m)*(n-m-1)/(2*(2*n-1))*(cos(theta))^(n-m-2) ...
                    + (n-m)*(n-m-1)*(n-m-3)/(8*(2*n-1)*(2*n-3))*(cos(theta))^(n-m-4));
            newTerm_inner = (g(n,mInd)*cos(m*lambda) + h(n,mInd)*sin(m*lambda)) * P;
            innerTerm = innerTerm + newTerm_inner;
        end

        newTerm_outer = (RE/R).^(n+1) * innerTerm; %Q: elementwise of nonelementwise vector?
        V_sum = V_sum + newTerm_outer;

    end

    V = RE * V_sum;

end