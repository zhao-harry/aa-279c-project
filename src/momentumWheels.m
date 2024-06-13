function [Mw,Lwdot] = momentumWheels(Lw,Mc,control,w)
    A = control.A;
    Lwdot = pinv(A) * (-Mc - cross(w,A * Lw));

    % Check saturation and apply correction
    if any(abs(Lw) > control.LMaxWheel)
        Lwdot = (abs(Lw) > control.LMaxWheel).* Lwdot;
    end

    % Recalculate Mc
    Mw = -A*Lwdot - cross(w,A * Lw);
end