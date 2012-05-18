function [phi]=vanLeer(r)
    % Gives the vanLeer limiter functional values
    %
    % [phi]=vanLeer(r)
    %
    % phi: the vanLeer limiter functional values
    % r: Sweby's r factor
    
    % Calculation
    phi=(r+abs(r))./(1+abs(r));
    
end


