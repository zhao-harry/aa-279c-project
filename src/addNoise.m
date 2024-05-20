function xNoise = addNoise(x, xRange)
    noiseDimensions = size(x);
    noise = rand(noiseDimensions);
    
    xNoise = x + (2 * xRange .* noise - xRange);
end