function A = DAD(m, v)
    % Implement the Deterministic Atittude Determination algorithm
    % Inputs:
    % - m: measured vectors ([m1 m2 . . . ])
    % - v: ground truth vector ([v1 v2 . . . ])
    % Outputs:
    % - A: DCM between principal axes and ECI

    A = m * pinv(v);

end