function A = DAD(m1, m2, v1, v2)
    % Implement the Deterministic Atittude Determination algorithm
    % Inputs:
    % - m1, m2: measured vectors
    % - v1, v2: ground truth vector
    % Outputs:
    % - A: DCM between principal axes and ECI

    m1_tilde = (m1 + m2)/2;
    m2_tilde = (m1 - m2)/2;

    pm = m1_tilde;
    qm = cross(m1_tilde,m2_tilde)/norm(cross(m1_tilde,m2_tilde));
    rm = cross(pm, qm);

    pv = v1;
    qv = cross(v1,v2)/norm(cross(v1,v2));
    rv = cross(pv,qv);

    M = [pm qm rm];
    V = [pv qv rv];
    A = M * pinv(V);
end