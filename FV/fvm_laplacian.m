function [A, RHS]=fvm_laplacian(field,gammaf,xC,xF,Sf)
    % Gives the matrix and corresponding RHS for div(gamma*grad(field))
    % 
    % [A, RHS]=fvm_laplacian(field,gamma,xC,XF)
    %
    % A: the matrix for discretized operator
    % RHS: right hand side column vector
    % field: the field the operator is applied to
    % gammaf: face diffusivity
    % xC: cell centres
    % xF: face centres
    % Sf: face areas

    % Recover field size
    N=size(field.internal,1);

    % Vectorized assembling of diffusion matrix
    % Diagonal part of advection matrix from 2 to N-1 (i.e. w/o BC)
    A=-diag([0; Sf(2:end-2).*gammaf(2:end-2)./(xC(2:end-1)-xC(1:end-2))+Sf(3:end-1).*gammaf(3:end-1)./(xC(3:end)-xC(2:end-1)); 0]);
      
    % Off diagonal parts of diffusion matrix
    A=A+diag([0; Sf(3:end-1).*gammaf(3:end-1)./(xC(3:end)-xC(2:end-1))],+1)+...
        diag([Sf(2:end-2).*gammaf(2:end-2)./(xC(2:end-1)-xC(1:end-2)); 0],-1);
    
    % Apply BC's
    RHS=zeros(N,1);
    
    % Left BC
    if (field.left.type=='V')
	  RHS(1)=-field.left.setvalue*gammaf(1)*Sf(1)/(xC(1)-xF(1));
	  A(1,1)=-(gammaf(1)*Sf(1)/(xC(1)-xF(1))+gammaf(2)*Sf(2)/(xC(2)-xC(1)));
	  A(1,2)=gammaf(2)*Sf(2)/(xC(2)-xC(1));
    elseif (field.left.type=='BP' | field.left.type=='G')
	  RHS(1)=gammaf(1)*field.left.gradient*Sf(1);
	  A(1,1)=-gammaf(2)*Sf(2)/(xC(2)-xC(1));
	  A(1,2)=gammaf(2)*Sf(2)/(xC(2)-xC(1));
    end
    
    % Right BC  
    if (field.right.type=='V')
	  RHS(N)=-field.right.setvalue*gammaf(end)*Sf(end)/(xF(end)-xC(end));
	  A(N,N)=-(gammaf(end)*Sf(end)/(xF(end)-xC(end))+gammaf(end-1)*Sf(end-1)/(xC(end)-xC(end-1)));
	  A(N,N-1)=gammaf(end-1)*Sf(end-1)/(xC(end)-xC(end-1));
    elseif (field.left.type=='BP' | field.left.type=='G')
	  RHS(end)=-gammaf(end)*field.right.gradient*Sf(end);
	  A(N,N)=-gammaf(end-1)*Sf(end-1)/(xC(end)-xC(end-1));
	  A(N,N-1)=gammaf(end-1)*Sf(end-1)/(xC(end)-xC(end-1));
    end
end
