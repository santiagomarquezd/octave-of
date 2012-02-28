function [Uf]=fvc_interpolate(U, w, xC, xF)
    % Gives a face linear interpolated field from cell centered field
    % 
    % [Uf]=fvm_interpolate(U, weights)
    %
    % U: cell centered field
    % w: linear interpolation weights
    % Uf: face interpolated field
    % xC: cell centers
    % xF: face centers

	% Face interpolated field allocation
	Uf=zeros(size(U.internal,1)+1,1);

	% Internal faces interpolation
	Uf(2:end-1)=w.*U.internal(1:end-1)+(1-w).*U.internal(2:end);	
	
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
