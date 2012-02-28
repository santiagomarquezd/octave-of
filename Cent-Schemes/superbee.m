function [phi]=superbee(r)
    % Gives the Superbee limiter functional values
    %
    % [phi]=superbee(r)
    %
    % phi: the Superbee limiter functional values
    % r: Sweby's r factor
    
    % Calculation
    phi=max(0,max(min(2*r,1),min(r,2)));
    
end


