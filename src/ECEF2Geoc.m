function geoc = ECEF2Geoc(rECEF,RE)
    % RE = 6378.1;
    lat = asind ( rECEF(3) / norm(rECEF) ) ;
    lon = atan2d ( rECEF(2) , rECEF(1) ) ;
    alt = norm ( rECEF ) - RE ;
    geoc = [ lat , lon , alt ];
end