function A = DAD(m1, m2, v1, v2)
    % Implement the Deterministic Atittude Determination algorithm
    % Inputs:
    % - m1, m2: measured vectors
    % - v1, v2: ground truth vector
    % Outputs:
    % - A: DCM between principal axes and ECI

    pm = m1;
    qm = cross(m1,m2)/norm(cross(m1,m2));
    rm = cross(pm, qm);

    pv = v1;
    qv = cross(v1,v2)/norm(cross(v1,v2));
    rv = cross(pv,qv);

    M = [pm qm rm];
    V = [pv qv rv];
    A = M * pinv(V);
end