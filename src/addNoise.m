function xNoise = addNoise(x, xRange, bias, biasRange)
    if nargin < 3
        bias = 0;
        biasRange = 0;
    end
    noiseDimensions = size(x);
    noise = randn(noiseDimensions);
    biasNoise = randn(noiseDimensions);
    xBias = bias*ones(noiseDimensions);
    
    xNoise = x + (2 * xRange .* noise - xRange) + xBias;
    % xNoise = x + (noise .* xRange) + (biasNoise .* biasRange) .* xBias;
end