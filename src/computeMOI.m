function I = computeMOI(filename,reference)
    data = readmatrix(filename);
    x = data(:,1) - reference(1);
    y = data(:,2) - reference(2);
    z = data(:,3) - reference(3);
    m = data(:,4);
    m_bus = m(1);
    m_RIS = m(2);
    m_panel = m(3);
    m_RAB = m(5);
    m_RAR = m(6);

    L_bus = 1.2;
    W_bus = 1.8;
    H_bus = 1.9;
    I_bus = m_bus * [(W_bus^2 + H_bus^2) / 12, 0, 0; ...
        0, (L_bus^2 + H_bus^2) / 12, 0; ...
        0, 0, (L_bus^2 + W_bus^2) / 12];

    L_RIS = 2.5;
    S_RIS = 0.459;
    a_RIS = S_RIS / (2 * tan(deg2rad(22.5)));
    R_avg = mean([a_RIS sqrt(a_RIS^2 + (S_RIS / 2)^2)]);
    I_RIS = m_RIS * [S_RIS^2 / 24 + a_RIS^2 / 2, 0, 0; ...
        0, L_RIS^2 / 12 + R_avg^2 / 4, 0; ...
        0, 0, L_RIS^2 / 12 + R_avg^2 / 4];

    W_panel = 6;
    H_panel = 1.9;
    I_panel = m_panel * [(W_panel^2 + H_panel^2) / 12, 0, 0; ...
        0, H_panel^2 / 12, 0; ...
        0, 0, W_panel^2 / 12];

    L_RAB = 0.1778;
    W_RAB = 0.1778;
    H_RAB = 9;
    deg_RAB = -18;
    rot_RAB = [cosd(deg_RAB), 0, sind(deg_RAB); ...
        0, 1, 0; ...
        -sind(deg_RAB), 0, cosd(deg_RAB)];
    I_RAB = m_RAB * [(W_RAB^2 + H_RAB^2) / 12, 0, 0; ...
        0, (L_RAB^2 + H_RAB^2) / 12, 0; ...
        0, 0, (L_RAB^2 + W_RAB^2) / 12];
    I_RAB_rot = rot_RAB * I_RAB * rot_RAB';

    R_RAR = 6;
    deg_RAR = -3.87;
    rot_RAR = [cosd(deg_RAR), 0, sind(deg_RAR); ...
        0, 1, 0; ...
        -sind(deg_RAR), 0, cosd(deg_RAR)];
    I_RAR = m_RAR * [R_RAR^2 / 4, 0, 0; ...
        0, R_RAR^2 / 4, 0; ...
        0, 0, R_RAR^2 / 2];
    I_RAR_rot = rot_RAR * I_RAR * rot_RAR';

    I_c = {I_bus, I_RIS, I_panel, I_panel, I_RAB_rot, I_RAR_rot};
    I = zeros([3 3]);
    for i = 1:length(I_c)
        D = [y(i)^2 + z(i)^2, -x(i) * y(i), -x(i) * z(i); ...
            -y(i) * x(i), x(i)^2 + z(i)^2, -y(i) * z(i); ...
            -z(i) * x(i), -y(i) * z(i), x(i)^2 + y(i)^2];
        I = I + I_c{i} + m(i) * D;
    end
end