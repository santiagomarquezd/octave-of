function [avalue]=Aop(M, V)
    % Gives the A matrix operator used in PISO loop
    % 
    % [avalue]=Aop(M, V)
    %
    % avalue: column vector of A opearator values
    % M: system matrix
    % V: column vector of volumes

    avalue=diag(M)./V;

end
