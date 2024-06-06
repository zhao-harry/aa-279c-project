function eulerAngs = A2eVec(A)
    len = size(A,3);
    eulerAngs = nan(3, len);
    for n = 1:len
        eulerAngs(:,n) = A2e(A(:,:,n));
    end
end