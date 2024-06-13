function [mode,m,Mmag] = magnetorquer(mode,B,Lw,control)
    % Check if we need to dump momentum
    if all(abs(Lw) < control.LwLimit * 0.1) && mode == 1
        mode = 0; % Stop desaturation
    end
    if any(abs(Lw) > control.LwLimit) || mode == 1
        mode = 1; % Trigger desaturation
        m = control.KMag / norm(B) * cross(control.A*Lw, normVec(B));
        
        % Implement magnetorquer saturation
        if any(m > control.magMax)
            magOver = m(m > control.magMax);
            scaleOver = max(magOver ./ control.magMax);
            m = m ./ scaleOver;
        end

        Mmag = cross(m,B);
    else
        mode = 0;
        m = [0;0;0];
        Mmag = [0;0;0];
    end
end
