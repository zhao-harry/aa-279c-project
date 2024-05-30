function measurementStarTracker = starTracker(qkminus,groundTruthVectors)
    % Extract information
    A_ECI2P = q2A(qkminus);
    rStarECI = groundTruthVectors(:,1:end-1);
    rStarP = A_ECI2P * rStarECI;

    measurementStarTracker = normVec(rStarP);
end