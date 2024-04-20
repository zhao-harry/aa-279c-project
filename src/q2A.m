function A = q2A(q)
    qv = q(1:3)';
    q4 = q(4);
    qx = [0, -qv(3), qv(2); ...
        qv(3), 0, -qv(1); ...
        -qv(2), qv(1), 0];
    A = (q4^2 - norm(qv)^2) * eye(3) + 2 * (qv * qv') - 2 * q4 * qx;
end