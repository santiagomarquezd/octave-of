function [phiR] = relaxF(phi, phi0, omega)
    % Relaxes field by omega factor
    %
    % [N, NRHS] = relaxM(M, RHS, omega)
    %
    % phi: update field value
    % phi0: previous step field value
    % omega: relaxation factor
    
    phiR = omega*phi.internal + (1 - omega)*phi0.internal;

end

