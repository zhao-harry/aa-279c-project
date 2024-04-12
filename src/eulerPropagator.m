function [wx, wy, wz] = eulerPropagator(Ix, Iy, Iz, tf, dt, wx, wy, wz)
    
    % Iteration
    for n = 1:(tf/dt)
        wx_dot = (Iy-Iz)/Ix*wy*wz;
        wy_dot = (Iz-Ix)/Iy*wz*wx;
        wz_dot = (Ix-Iy)/Iz*wx*wy;

        wx = wx_dot*dt + wx;
        wy = wy_dot*dt + wy;
        wz = wz_dot*dt + wz;
    end
end