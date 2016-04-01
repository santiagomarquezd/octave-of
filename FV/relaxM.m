function [N, NRHS] = relaxM(M, RHS, phi0, omega)
    % Relaxes matrix by omega factor
    %
    % [N, NRHS] = relaxM(M, RHS, omega)
    %
    % M: matrix
    % RHS: right hand side 
    % phi0: previous step field value
    % omega: relaxation factor
    
    N = M + diag(diag(M)*(1 - omega)/omega);
    NRHS = RHS + (1 - omega)/omega.*diag(M).*phi0.internal;

end

