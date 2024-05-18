function xNoise = addNoise(x, xRange)
    xNoise = x - xRange + 2*xRange*rand(size(x));
end