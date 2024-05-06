function GMST = time2GMST(time,MJD0)
    % Returns GMST in radians (wrapped to 2 pi)
    % MJD0 is MJD at start of simulation
    % Time is simulation time
    MJD = MJD0 + time / 86400;
    d = MJD - 51544.5;
    GMST = wrapTo2Pi(deg2rad(280.4606 + 360.9856473 * d));
end