function measurementSunSensor = sunSensor(qkminus,groundTruthVectors)
    % Get modeled parameters
    rSCI = groundTruthVectors(:,end);
    A_ECI2P = q2A(qkminus);

    % Get sensor measurement without noise
    rSunSensors = A_ECI2P * rSCI;
    measurementSunSensor = normVec(rSunSensors);
end