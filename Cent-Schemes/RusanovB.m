function [u1,u2]=RusanovB(u1_0,u2_0,flux1,flux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,dx,dt)
    % Applies the Rusanov scheme given actual state and cell fluxes
    %
    % [u1,u2]=Rusanov(u1_0,u2_0,flux1,flux2,cflux1,cflux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,dx,dt)
    %
    % u1,u2: new states for u1 and u2
    % u1_0,u2_0: previous states for u1 and u2
    % flux1,flux2: cell fluxes for u1 and u2 (with flux at boundaries as BC's)
    % cflux1,cflux2: conservative face fluxes from PISO loop
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % S1, S2: vector of source terms
    % rhom0: mixture density at previous timestep (needed when rhom factor is present in temporal derivative)
    % rhom: mixture density at actual timestep (needed when rhom factor is present in temporal derivative)
    % dx: spatial step
    % dt: time-step
    
    % Precomputing
    lambda=dt/dx;
    
    % Allocation
    u1=u1_0;
    u2=u2_0;
   
    % Rosunov scheme for non boundary cells (u1)
    u1.internal(2:end-1)=u1_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-lambda/2.*(flux1(3:end-1)-flux1(2:end-2))./rhom.internal(2:end-1)+1/2*lambda.*(a_j_plus_half(2:end-1).*(u1_0.internal(3:end)-u1_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)))./rhom.internal(2:end-1);

    % BC's
    u1.internal(1)=u1_0.internal(1).*rhom0.internal(1)./rhom.internal(1)-lambda*(flux1(2)-flux1(1))./rhom.internal(1)+lambda/2*(a_j_plus_half(1)*(u1_0.internal(2)-u1_0.internal(1))-a_j_minus_half(1)*(u1_0.internal(1)-u1_0.left.setvalue))./rhom.internal(1);
    u1.internal(end)=u1_0.internal(end).*rhom0.internal(end)./rhom.internal(end)-lambda*(flux1(end)-flux1(end-1))./rhom.internal(end)+lambda/2*(a_j_plus_half(end)*(u1_0.right.setvalue-u1_0.internal(end))-a_j_minus_half(end)*(u1_0.internal(end)-u1_0.internal(end-1)))./rhom.internal(end);
    
    u1.internal=u1.internal+S1*dt;
    u2.internal=u1.internal;

    % Putting second field is left


end
