function [phi]=minmod(r)
    % Gives the minmod limiter functional values
    %
    % [phi]=minmod(r)
    %
    % phi: the minmod limiter functional values
    % r: Sweby's r factor
    
    % Calculation
    phi=max(0,min(1,r));
    
end


