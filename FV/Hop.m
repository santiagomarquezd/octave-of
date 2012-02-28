function [hvalue]=Hop(M, RHS, field, V)
    % Gives the H matrix operator used in PISO loop
    % 
    % [hvalue]=Hop(M, RHS, V)
    %
    % hvalue: column vector of H opearator values
    % M: system matrix
    % RHS: right hand side of the system
    % field: solution field for M\RHS
    % V: column vector of volumes
    
    hvalue=(-(M-diag(diag(M)))*field.internal+RHS)./V;

end
