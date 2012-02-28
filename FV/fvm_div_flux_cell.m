function [A, RHS]=fvm_div_flux_cell(phi,field,xC,xF,weights,method)
    % Gives the matrix and corresponding RHS for divergence operator
    % from face flux and cell centered field
    % sum_f (phi)_f*(T)_f
    % 
    % fvm_div_flux_cell(phi,Sf,method)
    %
    % phi: flux at faces
    % field: the field T in div(U*T) calculation
    % xC: cell centers
    % xF: face centers
    % weights(2:end).*phi(3:end-1)
    % method: 1. Full upwind, 2. Linear interpolation

    % Recover field size
    N=size(field.internal,1);

    if (method==1)
      % Full upwind
      % Vectorized assembling of advection matrix
      % Diagonal part of advection matrix
      % A=diag(pos(phi(2:end)).*phi(2:end))+diag(neg(phi(1:end-1)).*phi(1:end-1));
      A=diag(pos(phi(2:end)).*phi(2:end))-diag(neg(phi(1:end-1)).*phi(1:end-1));
      
      % Off diagonal part of advection matrix
      % A=A-diag(pos(phi(2:end-1)).*phi(2:end-1),-1)-diag(neg(phi(2:end-1)).*phi(2:end-1),+1);
      % A=A+diag(pos(phi(2:end-1)).*phi(2:end-1),-1)+diag(neg(phi(2:end-1)).*phi(2:end-1),+1);
      A=A-diag(pos(phi(2:end-1)).*phi(2:end-1),-1)+diag(neg(phi(2:end-1)).*phi(2:end-1),+1);
            
      % Allocate RHS
      RHS=zeros(N,1);
      
      % Left BC
      if (field.left.type=='V')
		if (pos(phi(1)))
		  RHS(1,1)=phi(1)*field.left.setvalue;
		end
	  elseif (field.left.type=='G' | field.left.type=='BP')
		if (pos(phi(1)))
		  disp('WARNING!! Inconsistent left BC and divergence scheme (Fixed gradient in domain and Full upwind)')
		  % This kind of BC also affects the matrix
		  A(1,1)-=phi(1);
		  RHS(1,1)=phi(1)*field.left.gradient*abs(xC(1)-xF(1));
		elseif (neg(phi(1)))
		  % This kind of BC also affects the matrix
		  A(1,1)-=phi(1);
		  RHS(1,1)=phi(1)*field.left.gradient*abs(xC(1)-xF(1));
		end
	  end
		
	  % Right BC
      if (field.right.type=='V')
		if (neg(phi(end)))
		  RHS(end,1)=phi(end)*field.right.setvalue;
		end
	  elseif (field.right.type=='G' | field.right.type=='BP')
		if (neg(phi(end)))
		  disp('WARNING!! Inconsistent right BC and divergence scheme (Fixed gradient in domain enter and Full upwind)')
		  % This kind of BC also affects the matrix
		  A(end,end)+=phi(end);
		  RHS(end,1)=-phi(end)*field.right.gradient*abs(xF(end)-xC(end));
		elseif (pos(phi(end)))
		  % This kind of BC also affects the matrix
		  A(end,end)+=phi(end);
		  RHS(end,1)=-phi(end)*field.right.gradient*abs(xF(end)-xC(end));
		end
	  end	  

    elseif (method==2)
    
      % Linear interpolation
      % Vectorized assembling of advection matrix
      % Diagonal part of advection matrix
      A=diag([weights(1).*phi(2); weights(2:end).*phi(3:end-1)-(1-weights(1:end-1)).*phi(2:end-2); -(1-weights(end)).*phi(end-1)]);

      
      % Off diagonal part of advection matrix
      A=A+diag(-weights(1:end).*phi(2:end-1),-1)+diag((1-weights(1:end)).*phi(2:end-1),+1);

            
      % Allocate RHS
      RHS=zeros(N,1);
      
      % Left BC
      if (field.left.type=='V')
	RHS(1,1)=phi(1)*field.left.setvalue;
      elseif (field.left.type=='G' | field.left.type=='BP')
        % This kind of BC also affects the matrix
        A(1,1)-=phi(1);
	RHS(1,1)=phi(1)*field.left.gradient*abs(xC(1)-xF(1));
      end
		
      % Right BC
      if (field.right.type=='V')
	RHS(end,1)=-phi(end)*field.right.setvalue;
      elseif (field.left.type=='G' | field.left.type=='BP')
	% This kind of BC also affects the matrix
	A(end,end)+=phi(end);
	RHS(end,1)=-phi(end)*field.right.gradient*abs(xF(end)-xC(end));
      end	  
    end
end
