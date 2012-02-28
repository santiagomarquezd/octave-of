function [u1,u2]=Rusanov(u1_0,u2_0,flux1,flux2,a,dx,dt)
    % Applies the Rusanov scheme given actual state and fluxes
    %
    % [flux1,flux2]=no_Vr_cell_flux(U,alphag,rhom)
    %
    % u1,u2: new states for u1 and u2
    % u1_0,u2_0: previous states for u1 and u2
    % flux1,flux2: fluxes for u1 and u2
    % a: local speed of propagation at cell boundaries
    % dx: spatial step
    % dt: time-step
    
    % Precomputing
    lambda=dt/dx;
    
    % Allocation
    u1=u1_0;
    u2=u2_0;
    
    % Internal field caculation
    u1.internal(2:end-1)=u1_0.internal(2:end-1)-lambda./2.*(flux1.internal(3:end)-flux1.internal(1:end-2))+lambda./2.*(a(3:end-1).*(u1_0.internal(3:end)-u1_0.internal(2:end-1))-a(2:end-2).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)));
    u2.internal(2:end-1)=u2_0.internal(2:end-1)-lambda./2.*(flux2.internal(3:end)-flux2.internal(1:end-2))+lambda./2.*(a(3:end-1).*(u2_0.internal(3:end)-u2_0.internal(2:end-1))-a(2:end-2).*(u2_0.internal(2:end-1)-u2_0.internal(1:end-2)));
    
    % BC's
    u1.internal(1)=u1_0.internal(1)-lambda/2*(flux1.internal(2)-flux1.left.setvalue)+lambda/2*(a(2)*(u1_0.internal(2)-u1_0.internal(1))-a(1)*(u1_0.internal(1)-u1_0.left.setvalue));
    u1.internal(end)=u1_0.internal(end)-lambda/2*(flux1.right.setvalue-flux1.internal(end-1))+lambda/2*(a(end)*(u1_0.right.setvalue-u1_0.internal(end))-a(end-1)*(u1_0.internal(end)-u1_0.internal(end-1)));
    
    u2.internal(1)=u2_0.internal(1)-lambda/2*(flux2.internal(2)-flux2.left.setvalue)+lambda/2*(a(2)*(u2_0.internal(2)-u2_0.internal(1))-a(1)*(u2_0.internal(1)-u2_0.left.setvalue));
    u2.internal(end)=u2_0.internal(end)-lambda/2*(flux2.right.setvalue-flux2.internal(end-1))+lambda/2*(a(end)*(u2_0.right.setvalue-u2_0.internal(end))-a(end-1)*(u2_0.internal(end)-u2_0.internal(end-1)));
end

