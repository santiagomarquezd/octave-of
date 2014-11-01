function [k1,k2,k3]=RKLxFPoly(u1_0,u2_0,u3_0,V0,a,alphaDPL,dx,dt)
    % Applies the Lax-Friederichs scheme for polydisperse cases
    %
    % [u1,u2]=LxFPoly(u1_0,u2_0,flux1,flux2,dx,dt)
    %
    % u1,u2,u3: new states (only internal field)
    % u1_0,u2_0,u3_0: previous states (only internal field)
    % V0: relative velocity law constant
    % a: exponents vector
    % alphaDPL: concentration of the dense-packed layer
    % dx: spatial step
    % dt: time-step
 
    % Data reshaping
    u=[ (u1_0)' 
        (u2_0)' 
        (u3_0)'];

    F=arrayPFlux(u,V0,a,alphaDPL);
    % Face fluxes (impermeable walls)
    F=[[0; 0; 0] (F(:,1:(end-1))+F(:,2:(end)))/2 [0; 0; 0]];

    % LxF fluxes (the boundary fluxes are left zero)
    u=[u1_0 u2_0 u3_0]'; 
    F(:,2:end-1)=F(:,2:end-1)-1/2.*dx/dt.*(u(:,2:end)-u(:,1:end-1));

    % Rusanov integration
    k=-1/dx*(F(:,2:end)-F(:,1:end-1));

    k1=k(1,:)';
    k2=k(2,:)';
    k3=k(3,:)';

end

