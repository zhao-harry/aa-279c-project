function eulerAngs = quats2Euler(quats)
    nLen = size(quats, 2);
    eulerAngs = nan(3,nLen);
    for n = 1:nLen
        eulerAngs(:,n) = A2e(q2A(quats(:,n)));
    end
end