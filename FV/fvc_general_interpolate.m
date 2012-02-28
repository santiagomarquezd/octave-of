function [Uf]=fvc_general_interpolate(U, xC, xF, scheme, directionFlux)
    % Gives a general face interpolated field from cell centered field
    % 
    % [Uf]=fvm_general_interpolate(U, xC, xF, scheme, directionFlux)
    %
    % U: cell centered field
    % Uf: face interpolated field
    % xC: cell centers
    % xF: face centers
    % scheme: -1. Upwind, 0. Linear. 1, Downwind
    % directionFlux: flux giving the flow direction (+1, -1 in each
    % INTERNAL face)
    %
    
	% Face interpolated field allocation
	Uf=zeros(size(U.internal,1)+1,1);
	
	if (scheme==-1)
	  % Internal faces interpolation
	  Uf(2:end-1)=U.internal(1:end-1).*(directionFlux+1)./2+...
	  +U.internal(2:end).*(directionFlux-1)./(-2);	
	elseif (scheme==0)
	      % Interpolation weights calculation
	      w=weights(xC, xF);
	      Uf(2:end-1)=w.*U.internal(1:end-1)+(1-w).*U.internal(2:end);
	elseif (scheme==1)
	  % Internal faces interpolation
	  Uf(2:end-1)=U.internal(1:end-1).*(directionFlux-1)./(-2)+...
	  +U.internal(2:end).*(directionFlux+1)./(2);	
	end
	
	% Boundary conditions
	if (U.left.type=='V')
	  Uf(1)=U.left.value;
	elseif (U.left.type=='G')
	  Uf(1)=U.internal(1)-U.left.gradient*(xC(1)-xF(1));
	end

	if (U.right.type=='V')
	  Uf(end)=U.right.value;
	elseif (U.right.type=='G')
	  Uf(end)=U.internal(end)+U.right.gradient*(xF(end)-xC(end));
	end


end
