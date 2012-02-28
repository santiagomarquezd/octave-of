function [phiC] = fvc_ddtPhiCorrection(U0, phi0, rUA, rUAf, rho0, w, xC, xF, Sf, dt)
    % Explicit correction of mass flux taken into account changes in time
    %
    % [phiC] = fvc_ddtPhiCorrection(U0, phi0, rUA, rUAf, rho0, w, xC, xF, Sf, dt)
    %
    % U0: velocity at previous timestep
    % phi0: volumetric flux at previous timestep
    % rUA: 1/Aoperator
    % rUAf: rUA at faces
    % rho0: density at previous timestep
    % w: interpolation weights
    % xC: cell centres positions
    % xF: face centres positions
    % Sf: face areas
    % dt: timestep  

    Kc=fvc_ddtPhiCoeff(U0, phi0./fvc_interpolate(rho0,w,xC,xF), w, xC, xF, Sf);

    phiC=1/dt.*Kc.*(fvc_interpolate(assign(rUA,rho0,'*'),w,xC,xF).* ...
	  phi0./fvc_interpolate(rho0,w,xC,xF)-...
	  (fvc_interpolate(assign(assign(rUA,rho0,'*'),U0,'*'),w,xC,xF).*Sf));
	  
end





