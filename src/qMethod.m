function q = qMethod(m, v, w)
    % Implement the Deterministic q method
    % Inputs:
    % - m: measured vectors
    % - v: ground truth vector
    % - w: weights of sensor measurements
    % Outputs:
    % - A: DCM between principal axes and ECI

    U = sqrt(w) .* v;
    W = sqrt(w) .* m;

    B = W * U';
    S = B + B';
    Z = [B(2,3)-B(3,2), B(3,1)-B(1,3), B(1,2)-B(2,1)]';
    sigma = trace(B);

    K = [S - eye(size(S))*sigma, Z;
            Z', sigma];

    [V, lamda] = eig(K);

    [~, nMax] = max(diag(lamda));

    q = V(:, nMax);
    q = q * q(1)/abs(q(1));
end