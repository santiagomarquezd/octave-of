function [u1,u2]=Rusanov(u1_0,u2_0,flux1,flux2,cflux1,cflux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,dx,dt)
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
   
    if 0

      	  % Rosunov scheme for non boundary cells (u1)
	  u1.internal(2:end-1)=u1_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-lambda/2./rhom.internal(2:end-1).*(flux1.internal(3:end)-flux1.internal(1:end-2))+1/2*lambda.*(a_j_plus_half(2:end-1).*(u1_0.internal(3:end)-u1_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)))-1./rhom.internal(2:end-1).*(cflux1(3:end-1)-cflux1(2:end-2))*lambda;
   
	  % BC's
	  u1.internal(1)=u1_0.internal(1)*rhom0.internal(1)/rhom.internal(1)-lambda/rhom.internal(1)*(flux1.internal(2)-flux1.left.setvalue)+lambda/2*(a_j_plus_half(1)*(u1_0.internal(2)-u1_0.internal(1)))-lambda/rhom.internal(1)*(a_j_minus_half(1)*(u1_0.internal(1)-u1_0.left.setvalue))-1./rhom.internal(1)*(cflux1(2)-cflux1(1))*lambda;
	  u1.internal(end)=u1_0.internal(end)*rhom0.internal(end)/rhom.internal(end)-lambda/rhom.internal(end)*(flux1.right.setvalue-flux1.internal(end-1))+lambda*(a_j_plus_half(end)*(u1_0.right.setvalue-u1_0.internal(end)))-lambda/rhom.internal(end)/2*(a_j_minus_half(end)*(u1_0.internal(end)-u1_0.internal(end-1)))-1./rhom.internal(end)*(cflux1(end)-cflux1(end-1))*lambda;
	  
	  u1.internal=u1.internal+S1*dt./rhom.internal;

	  %keyboard; pause;

	  % Rosunov scheme for non boundary cells (u2)
	  u2.internal(2:end-1)=u2_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-lambda/2./rhom.internal(2:end-1).*(flux2.internal(3:end)-flux2.internal(1:end-2))+1/2*lambda.*(a_j_plus_half(2:end-1).*(u2_0.internal(3:end)-u2_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u2_0.internal(2:end-1)-u2_0.internal(1:end-2)))-1./rhom.internal(2:end-1).*(cflux2(3:end-1)-cflux2(2:end-2))*lambda;
   
	  % BC's
	  u2.internal(1)=u2_0.internal(1)*rhom0.internal(1)/rhom.internal(1)-lambda/rhom.internal(1)*(flux2.internal(2)-flux2.left.setvalue)+lambda/2*(a_j_plus_half(1)*(u2_0.internal(2)-u2_0.internal(1)))-lambda/rhom.internal(1)*(a_j_minus_half(1)*(u2_0.internal(1)-u2_0.left.setvalue))-1./rhom.internal(1)*(cflux2(2)-cflux2(1))*lambda;
	  u2.internal(end)=u2_0.internal(end)*rhom0.internal(end)/rhom.internal(end)-lambda/rhom.internal(end)*(flux2.right.setvalue-flux2.internal(end-1))+lambda*(a_j_plus_half(end)*(u2_0.right.setvalue-u2_0.internal(end)))-lambda/2/rhom.internal(end)*(a_j_minus_half(end)*(u2_0.internal(end)-u2_0.internal(end-1)))-1./rhom.internal(end)*(cflux2(end)-cflux2(end-1))*lambda;
 
          u2.internal=u2.internal+S2*dt./rhom.internal;
    
    else

	  %keyboard; pause;

	  % Rosunov scheme for non boundary cells (u1)
	  u1.internal(2:end-1)=u1_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-...
			       lambda/2./rhom.internal(2:end-1).*(flux1.internal(3:end)-flux1.internal(1:end-2))+...
			       1/2*lambda./rhom.internal(2:end-1).*rhom0.internal(2:end-1).*(a_j_plus_half(2:end-1).*(u1_0.internal(3:end)-u1_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)))-...
			       1./rhom.internal(2:end-1).*(cflux1(3:end-1)-cflux1(2:end-2))*lambda;
   
	  % BC's
	  u1.internal(1)=u1_0.internal(1)*rhom0.internal(1)/rhom.internal(1)-lambda/rhom.internal(1)*(flux1.internal(2)-flux1.left.setvalue)+lambda/2/rhom.internal(1)*rhom0.internal(1)*(a_j_plus_half(1)*(u1_0.internal(2)-u1_0.internal(1)))-lambda/rhom.internal(1)*(a_j_minus_half(1)*(u1_0.internal(1)-u1_0.left.setvalue))-1./rhom.internal(1)*(cflux1(2)-cflux1(1))*lambda;
	  u1.internal(end)=u1_0.internal(end)*rhom0.internal(end)/rhom.internal(end)-lambda/rhom.internal(end)*(flux1.right.setvalue-flux1.internal(end-1))+lambda/rhom.internal(end)**rhom0.internal(end)*(a_j_plus_half(end)*(u1_0.right.setvalue-u1_0.internal(end)))-lambda/rhom.internal(end)/2*(a_j_minus_half(end)*(u1_0.internal(end)-u1_0.internal(end-1)))-1./rhom.internal(end)*(cflux1(end)-cflux1(end-1))*lambda;
	  
	  u1.internal=u1.internal+S1*dt./rhom.internal;

	  %keyboard; pause;

	  % Rosunov scheme for non boundary cells (u2)
	  u2.internal(2:end-1)=u2_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-...
			       lambda/2./rhom.internal(2:end-1).*(flux2.internal(3:end)-flux2.internal(1:end-2))+...
			       1/2*lambda./rhom.internal(2:end-1).*rhom0.internal(2:end-1).*(a_j_plus_half(2:end-1).*(u2_0.internal(3:end)-u2_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u2_0.internal(2:end-1)-u2_0.internal(1:end-2)))-...
			       1./rhom.internal(2:end-1).*(cflux2(3:end-1)-cflux2(2:end-2))*lambda;
   
	  % BC's
	  u2.internal(1)=u2_0.internal(1)*rhom0.internal(1)/rhom.internal(1)-lambda/rhom.internal(1)*(flux2.internal(2)-flux2.left.setvalue)+lambda/2/rhom.internal(1)*rhom0.internal(1)*(a_j_plus_half(1)*(u2_0.internal(2)-u2_0.internal(1)))-lambda/rhom.internal(1)*(a_j_minus_half(1)*(u2_0.internal(1)-u2_0.left.setvalue))-1./rhom.internal(1)*(cflux2(2)-cflux2(1))*lambda;
	  u2.internal(end)=u2_0.internal(end)*rhom0.internal(end)/rhom.internal(end)-lambda/rhom.internal(end)*(flux2.right.setvalue-flux2.internal(end-1))+lambda/rhom.internal(end)*rhom0.internal(end)*(a_j_plus_half(end)*(u2_0.right.setvalue-u2_0.internal(end)))-lambda/2/rhom.internal(end)*(a_j_minus_half(end)*(u2_0.internal(end)-u2_0.internal(end-1)))-1./rhom.internal(end)*(cflux2(end)-cflux2(end-1))*lambda;
 
          u2.internal=u2.internal+S2*dt./rhom.internal;
    end
end

 %Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx*(phiAlpha(2:end).*...
  % 			 Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1))./rhom.internal(1:end);


  
