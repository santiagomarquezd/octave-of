function [u1,u2]=LxF(u1_0,u2_0,flux1,flux2,dx,dt)
    % Applies the Lax-Friederichs scheme given actual state and fluxes
    %
    % [flux1,flux2]=no_Vr_cell_flux(U,alphag,rhom)
    %
    % u1,u2: new states for u1 and u2
    % u1_0,u2_0: previous states for u1 and u2
    % flux1,flux2: fluxes for u1 and u2
    % dx: spatial step
    % dt: time-step
    
    % Allocation
    u1=u1_0;
    u2=u2_0;

    %keyboard; pause;
    
    % Internal field caculation
    u1.internal(2:end-1)=(u1_0.internal(1:end-2)+u1_0.internal(3:end))/2-dt/dx/2*(flux1.internal(3:end)-flux1.internal(1:end-2));
    u2.internal(2:end-1)=(u2_0.internal(1:end-2)+u2_0.internal(3:end))/2-dt/dx/2*(flux2.internal(3:end)-flux2.internal(1:end-2));
    
    % BC's
    u1.internal(1)=(u1_0.left.setvalue+u1_0.internal(2))/2-dt/dx/2*(flux1.internal(2)-flux1.left.setvalue);
    u1.internal(end)=(u1_0.right.setvalue+u1_0.internal(end-1))/2-dt/dx/2*(flux1.right.setvalue-flux1.internal(end-1));
    
    
    u2.internal(1)=(u2_0.left.setvalue+u2_0.internal(2))/2-dt/dx/2*(flux2.internal(2)-flux2.left.setvalue);
    u2.internal(end)=(u2_0.right.setvalue+u2_0.internal(end-1))/2-dt/dx/2*(flux2.right.setvalue-flux2.internal(end-1));
end