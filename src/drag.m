function [F,M] = drag(v,rho,CD,barycenter,normal,area,cm)
    % Compute drag in principal axes
    % RE = 6378.1 % km
    u = v / norm(v);
    N = normal;
    Aeff = ((u' * N) > 0) .* area;
    rC = barycenter - cm;
    D = -0.5 * CD * rho * norm(v)^2 * (u' * (N .* Aeff)) .* u;
    F = sum(D,2);
    M = sum(cross(rC,D),2);

    % Slow loop function (obsolete)
    % u = v / norm(v);
    % F = zeros([3 1]);
    % M = zeros([3 1]);
    % for i = length(area)
    %     n = normal(:,i);
    %     if dot(u,n) < 0
    %         rC = (barycenter(:,i) - cm);
    %         A = area(1,i);
    %         D = -0.5 * CD * rho * norm(v)^2 * dot(u,n) * u * A;
    %         M = M + cross(rC,D);
    %         F = F + D;
    %     end
    % end
end