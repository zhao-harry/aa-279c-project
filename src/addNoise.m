function xNoise = addNoise(x, xRange, bias)
    if nargin < 3
        bias = 0;
    end
    noiseDimensions = size(x);
    noise = randn(noiseDimensions);
    xBias = bias*ones(noiseDimensions);
    
    xNoise = x + (2 * xRange .* noise - xRange) + xBias;
end