function xNoise = addNoise(x, xRange)
    noiseDimensions = size(x);
    noise = normrnd(0, 0.5, noiseDimensions);
    
    xNoise = x + (2 * xRange .* noise - xRange);
end