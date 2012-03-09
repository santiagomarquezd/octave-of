function [u1,u2]=Rusanov(u1_0,u2_0,flux1,flux2,a_j_minus_half,a_j_plus_half,S1,S2,dx,dt)
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
   
	if 0
	  % Internal field caculation
	  u1.internal(2:end-1)=u1_0.internal(2:end-1)-lambda./2.*(flux1.internal(3:end)-flux1.internal(1:end-2))+lambda./2.*(a(3:end-1).*(u1_0.internal(3:end)-u1_0.internal(2:end-1))-a(2:end-2).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)));
	  u2.internal(2:end-1)=u2_0.internal(2:end-1)-lambda./2.*(flux2.internal(3:end)-flux2.internal(1:end-2))+lambda./2.*(a(3:end-1).*(u2_0.internal(3:end)-u2_0.internal(2:end-1))-a(2:end-2).*(u2_0.internal(2:end-1)-u2_0.internal(1:end-2)));
	  
	  % BC's
	  u1.internal(1)=u1_0.internal(1)-lambda*(flux1.internal(2)-flux1.left.setvalue)+lambda*(a(2)*(u1_0.internal(2)-u1_0.internal(1))-a(1)*(u1_0.internal(1)-u1_0.left.setvalue));
	  u1.internal(end)=u1_0.internal(end)-lambda*(flux1.right.setvalue-flux1.internal(end-1))+lambda*(a(end)*(u1_0.right.setvalue-u1_0.internal(end))-a(end-1)*(u1_0.internal(end)-u1_0.internal(end-1)));
	  
	  u2.internal(1)=u2_0.internal(1)-lambda*(flux2.internal(2)-flux2.left.setvalue)+lambda*(a(2)*(u2_0.internal(2)-u2_0.internal(1))-a(1)*(u2_0.internal(1)-u2_0.left.setvalue));
	  u2.internal(end)=u2_0.internal(end)-lambda*(flux2.right.setvalue-flux2.internal(end-1))+lambda*(a(end)*(u2_0.right.setvalue-u2_0.internal(end))-a(end-1)*(u2_0.internal(end)-u2_0.internal(end-1)));
    
    else
	  % Rosunov scheme for non boundary cells (u1)
	  u1(2:end-1)=u1_0(2:end-1)-lambda/2*(flux1(3:end)-flux1(1:end-2))+1/2*lambda*(a_j_plus_half(2:end-1).*(u1_0(3:end)-u1_0(2:end-1))-a_j_minus_half(2:end-1).*(u1_0(2:end-1)-u1_0(1:end-2)));
   
	  % BC's
	  u1(1)=u1_0(1)-lambda*(flux1(2)-flux1.left.setvalue)+lambda/2*(a_j_plus_half(1)*(u1_0(2)-u1_0(1)))-lambda*(a_j_minus_half(1)*(u1_0(1)-u1_0.left.setvalue));
	  u1(end)=u1_0(end)-lambda*(flux1.right.setvalue-flux1(end-1))+lambda*(a_j_plus_half(end)*(u1_0.right.setvalue-u1_0(end)))-lambda/2*(a_j_minus_half(end)*(u1_0(end)-u1_0(end-1)));
	  
	  % Rosunov scheme for non boundary cells (u2)
	  u2(2:end-1)=u2_0(2:end-1)-lambda/2*(flux2(3:end)-flux2(1:end-2))+1/2*lambda*(a_j_plus_half(2:end-1).*(u2_0(3:end)-u2_0(2:end-1))-a_j_minus_half(2:end-1).*(u2_0(2:end-1)-u2_0(1:end-2)));
   
	  % BC's
	  u2(1)=u2_0(1)-lambda*(flux2(2)-flux2.left.setvalue)+lambda/2*(a_j_plus_half(1)*(u2_0(2)-u2_0(1)))-lambda*(a_j_minus_half(1)*(u2_0(1)-u2_0.left.setvalue));
	  u2(end)=u2_0(end)-lambda*(flux2.right.setvalue-flux2(end-1))+lambda*(a_j_plus_half(end)*(u2_0.right.setvalue-u2_0(end)))-lambda/2*(a_j_minus_half(end)*(u2_0(end)-u2_0(end-1)));
 
    
    end
end


  
