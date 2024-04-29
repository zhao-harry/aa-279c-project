function w = eulerEquationWithWheelStepper(w0,Ix,Iy,Iz,Ir,tFinal,tStep)
    
    % Initialize
    nStep = ceil(tFinal/tStep);
    w = nan(nStep+1,4);
    w(1,:) = w0';
    wxDot = 0;
    wyDot = 0;
    
    for n = 1:nStep
        % Set parameters
        wx = w(n,1);
        wy = w(n,2);
        wz = w(n,3);
        wr = w(n,4);
        
        % Find derivatives
        a = (Iz - Iy)*wz + Ir*wr;
        b = (Ix - Iz)*wz - Ir*wr;
        wxDot = -a/Ix*dwy;
        wyDot = -b/Iy*dwx;

        % Propogate
        w(n+1,1) = wx + wxDot*tStep;
        w(n+1,2) = wy + wyDot*tStep;
        w(n+1,3) = wz;
        w(n+1,4) = wr;
    end
end