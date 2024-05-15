function A = DAD_twoVecs(m1, m2, v1, v2)
    % Implement the Deterministic Atittude Determination algorithm
    % Inputs:
    % - m1, m2: measured vectors
    % - v1, v2: ground truth vector
    % Outputs:
    % - A: DCM between principal axes and ECI

    m1_tilde = (m1 + m2)/2;
    m2_tilde = (m1 - m2)/2;
    v1_tilde = (v1 + v2)/2;
    v2_tilde = (v1 - v2)/2;

    m1_tilde = m1_tilde / norm(m1_tilde);
    m2_tilde = m2_tilde / norm(m2_tilde);
    v1_tilde = v1_tilde / norm(v1_tilde);
    v2_tilde = v2_tilde / norm(v2_tilde);
    
    pm = m1_tilde;
    cross_qm = cross(m1_tilde, m2_tilde);
    qm = cross_qm/norm(cross_qm);
    rm = cross(pm, qm);

    pv = v1;
    cross_qv = cross(v1_tilde, v2_tilde);
    qv = cross_qv/norm(cross_qv);
    rv = cross(pv,qv);

    M = [pm qm rm];
    V = [pv qv rv];
    A = M / V;
end