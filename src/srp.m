function [F,M] = srp(s,P,Cd,Cs,barycenter,normal,area,cm)
    % Compute solar radiation pressure in principal axes
    u = s / norm(s);
    N = normal;
    Aeff = ((u' * N) > 0) .* area;
    rC = barycenter - cm;
    theta = acos((u' * N) ./ (norm(u) * vecnorm(N)));
    SRP = -P * cos(theta) .* Aeff .* ...
        ((1 - Cs) * u + 2 * (Cs * cos(theta) + Cd / 3) .* N);
    F = sum(SRP,2);
    M = sum(cross(rC,SRP),2);
    
    % Slow loop function (obsolete)
    % F = zeros([3 1]);
    % M = zeros([3 1]);
    % for i = length(area)
    %     n = normal(:,i);
    %     if dot(u,n) > 0
    %         rC = barycenter(:,i) - cm;
    %         theta = acos(dot(u,n) / (norm(u) * norm(n)));
    %         A = area(1,i);
    %         SRP = -P * cos(theta) * A * ...
    %             ((1 - Cs) * u + 2 * (Cs * cos(theta) + Cd / 3) * n);
    %         M = M + cross(rC,SRP);
    %         F = F + SRP;
    %     end
    % end
end