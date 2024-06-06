function xNoise = addNoise(x, xRange, bias, biasRange)
    if nargin < 3
        bias = 0;
    end
    if nargin < 4
        biasRange = 0;
    end
    noiseDimensions = size(x);
    noise = randn(noiseDimensions);
    biasNoise = randn(noiseDimensions);
    xBias = bias*ones(noiseDimensions);
    
    % xNoise = x + (xRange .* noise) + xBias;
    xNoise = x + (noise .* xRange) + (biasNoise .* biasRange) .* xBias;
end