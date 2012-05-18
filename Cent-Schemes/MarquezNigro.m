function [u1]=MarquezNigro(u1_0,u2_0,flux1,flux2,cflux1,cflux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,w,xC,xF,dx,dt)
    % Applies the Marquez-Nigro scheme given actual state and cell fluxes
    %
    % [u1]=MarquezNigro(u1_0,u2_0,flux1,flux2,cflux1,cflux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,w,xC,xF,dx,dt)
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
    % w: linear interpolation weights
    % xC: cell centroids positions
    % xF: face centroids positions
    % dx: spatial step
    % dt: time-step
    
    % Precomputing
    lambda=dt/dx;
    
    if 1
      Flux1=cflux1.*fvc_interpolate(u1_0, w, xC, xF);
      Flux2=cflux2.*fvc_interpolate(u2_0, w, xC, xF);
    else
      directionFlux=sign(cflux1(2:end-1,1)./(abs(cflux1(2:end-1,1))+1E-9));
      Flux1=cflux1.*fvc_general_interpolate(u1_0, xC, xF,-1,directionFlux);
    end

    if 0
	  % Testing for boundary values using ghost cells (Leveque)
	  % The values of ghosts cell depend in the kind of BC

	  % Left BC
	  if (u1_0.left.type=='V')
	    % Fixed value, left ghost cell is the same of left BC
	    u1_0.internal=[u1_0.left.setvalue;u1_0.internal];
	    % For the conserved flux left ghost face flux is similar the real left face flux
	    cflux1=[cflux1(1);cflux1];
	    % The same for flux1 and Flux1
	    flux1.internal=[flux1.left.setvalue;flux1.internal];
	    Flux1=[Flux1(1);Flux1];
	  elseif (u1_0.left.type=='G')
	    % Works only for null gradient
	    % Fixed gradient, left ghost cell is the same if first interior cell
	    u1_0.internal=[u1_0.internal(1);u1_0.internal];
	    % For the conserved flux left ghost face flux is similar the real left face flux
	    cflux1=[cflux1(1);cflux1];
	    % The same for flux1 and Flux1
	    flux1.internal=[flux1.internal(1);flux1.internal];
	    Flux1=[Flux1(1);Flux1];
	  end

	  % Right BC
	  if (u1_0.right.type=='V')
	    % Fixed value, right ghost cell is the same of right BC
	    u1_0.internal=[u1_0.internal;u1_0.right.setvalue];
	    % For the conserved flux right ghost face flux is similar the real right face flux
	    cflux1=[cflux1;cflux1(end)];
	    % The same for flux1 and Flux1	      
	    flux1.internal=[flux1.internal;flux1.right.setvalue];
	    Flux1=[Flux1;Flux1(end)];
	  elseif (u1_0.right.type=='G')
	    % Works only for null gradient
	    % Fixed gradient, right ghost cell is the same if last interior cell
	    u1_0.internal=[u1_0.internal;u1_0.internal(end)];
	    % For the conserved flux right ghost face flux is similar the real right face flux
	    cflux1=[cflux1;cflux1(end)];
	    % The same for flux1 and Flux1	      
	    flux1.internal=[flux1.internal;flux1.internal(end)];
	    Flux1=[Flux1;Flux1(end)];
	  end

	  % In rhom and the speeds the setvalues field already has the correct values depending on the BC's (due the setBC)
	  rhom0.internal=[rhom0.left.setvalue;rhom0.internal;rhom0.right.setvalue];
	  rhom.internal=[rhom.left.setvalue;rhom.internal;rhom.right.setvalue];
	  a_j_minus_half=[a_j_minus_half(1);a_j_minus_half;a_j_plus_half(end)];
	  a_j_plus_half=[a_j_minus_half(1);a_j_plus_half;a_j_plus_half(end)];

	  % Allocation
	  u1=u1_0;

	  %keyboard; pause;

	  % Now, internal field and BC's are treated in the same calculus
	  u1.internal(2:end-1)=u1_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-...
			       lambda/2./rhom.internal(2:end-1).*(flux1.internal(3:end)-flux1.internal(1:end-2))+...
			       1/2*lambda./rhom.internal(2:end-1).*rhom0.internal(2:end-1).*(a_j_plus_half(2:end-1).*(u1_0.internal(3:end)-...
			       u1_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)))-...
			       1./rhom.internal(2:end-1).*(Flux1(3:end-1)-Flux1(2:end-2))*lambda+...
			       lambda^2/2./rhom.internal(2:end-1).*(a_j_plus_half(2:end-1).*((cflux1(3:end-1)-...
			       cflux1(2:end-2)).*u1_0.internal(2:end-1)-(cflux1(4:end)-...
			       cflux1(3:end-1)).*u1_0.internal(3:end))+a_j_minus_half(2:end-1).*((cflux1(3:end-1)-...
			       cflux1(2:end-2)).*u1_0.internal(2:end-1)-(cflux1(2:end-2)-...
			       cflux1(1:end-3)).*u1_0.internal(1:end-2)));
   	  
	  u1.internal=u1.internal+S1*dt./rhom.internal;

	  % Erasing ghost cells
	  u1.internal=u1.internal(2:end-1);
    
    else
	  % Allocation
	  u1=u1_0;
	  u2=u2_0;

	  %keyboard; pause;
	  u1.internal(2:end-1)=u1_0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-...
			       lambda/2./rhom.internal(2:end-1).*(flux1.internal(3:end)-flux1.internal(1:end-2))+...
			       1/2*lambda./rhom.internal(2:end-1).*rhom0.internal(2:end-1).*(a_j_plus_half(2:end-1).*(u1_0.internal(3:end)-...
			       u1_0.internal(2:end-1))-a_j_minus_half(2:end-1).*(u1_0.internal(2:end-1)-u1_0.internal(1:end-2)))-...
			       1./rhom.internal(2:end-1).*(Flux1(3:end-1)-Flux1(2:end-2))*lambda+...
			       lambda^2/2./rhom.internal(2:end-1).*(a_j_plus_half(2:end-1).*((cflux1(3:end-1)-...
			       cflux1(2:end-2)).*u1_0.internal(2:end-1)-(cflux1(4:end)-...
			       cflux1(3:end-1)).*u1_0.internal(3:end))+a_j_minus_half(2:end-1).*((cflux1(3:end-1)-...
			       cflux1(2:end-2)).*u1_0.internal(2:end-1)-(cflux1(2:end-2)-...
			       cflux1(1:end-3)).*u1_0.internal(1:end-2)));
   
	  % BC's
	  u1.internal(1)=u1_0.internal(1)*rhom0.internal(1)/rhom.internal(1)-...
			  lambda/rhom.internal(1)*(flux1.internal(2)-flux1.left.setvalue)+...
			  lambda/2/rhom.internal(1)*rhom0.internal(1)*(a_j_plus_half(1)*(u1_0.internal(2)-u1_0.internal(1)))-...
			  lambda/rhom.internal(1)*(a_j_minus_half(1)*(u1_0.internal(1)-u1_0.left.setvalue))-1./rhom.internal(1)*(Flux1(2)-...
			  Flux1(1))*lambda+...
			  2*lambda^2./rhom.internal(1).*(a_j_plus_half(1).*((cflux1(2)-...
			  cflux1(1)).*u1_0.internal(1)-(0))+a_j_minus_half(1).*((cflux1(2)-...
			  cflux1(1)).*u1_0.internal(1)-(0)));

	  u1.internal(end)=u1_0.internal(end)*rhom0.internal(end)/rhom.internal(end)-...
			  lambda/rhom.internal(end)*(flux1.right.setvalue-flux1.internal(end-1))+...
			  lambda/rhom.internal(end)*rhom0.internal(end)*(a_j_plus_half(end)*(u1_0.right.setvalue-u1_0.internal(end)))-...
			  lambda/rhom.internal(end)/2*(a_j_minus_half(end)*(u1_0.internal(end)-u1_0.internal(end-1)))-1./rhom.internal(end)*(Flux1(end)-...
			  Flux1(end-1))*lambda+...
			  2*lambda^2./rhom.internal(end).*(a_j_plus_half(end).*((cflux1(end)-...
			  cflux1(end-1)).*u1_0.internal(end)-(0))+a_j_minus_half(end).*((cflux1(end)-...
			  cflux1(end-1)).*u1_0.internal(end)-(0)));
	  
	  u1.internal=u1.internal+S1*dt./rhom.internal;

    end
end

% En linea de comando
% Alphag0.internal(2:end-1).*rhom0.internal(2:end-1)./rhom.internal(2:end-1)-lambda/2./rhom.internal(2:end-1).*(fluxAlpha.internal(3:end)-fluxAlpha.internal(1:end-2))+1/2*lambda./rhom.internal(2:end-1).*rhom0.internal(2:end-1).*(a_j_plus_half(2:end-1).*(Alphag0.internal(3:end)-Alphag0.internal(2:end-1))-a_j_minus_half(2:end-1).*(Alphag0.internal(2:end-1)-Alphag0.internal(1:end-2)))
% -1./rhom.internal(2:end-1).*(Flux1(3:end-1)-Flux1(2:end-2))*lambda+lambda^2/2./rhom.internal(2:end-1).*(a_j_plus_half(2:end-1).*((cFluxAlpha(3:end-1)-cFluxAlpha(2:end-2)).*Alphag0.internal(2:end-1)-(cFluxAlpha(4:end)-cFluxAlpha(3:end-1)).*Alphag0.internal(3:end))+a_j_minus_half(2:end-1).*((cFluxAlpha(3:end-1)-cFluxAlpha(2:end-2)).*Alphag0.internal(2:end-1)-(cFluxAlpha(2:end-2)-cFluxAlpha(1:end-3)).*Alphag0.internal(1:end-2)));  
