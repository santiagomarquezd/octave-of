function [out]=flux(M,RHS,field)
  % Gives the flux from matrix, RHS, and field
  %
  %
  % [out]=flux(M,RHS,field)
  % 
  % out: resulting flux
  % M: diffusive matrix
  % RHS: right hand side of the system
  % field: actual values of the field
  
  % out allocation
  out = zeros(size(field.internal,1)+1,1);
  
  % Internal faces flux calculation
  out(2:end-1)=diag(M,+1).*field.internal(2:end)-diag(M,-1).*field.internal(1:end-1);
  
  % Boundary faces
  if (field.left.type=='G' | field.left.type=='BP')
    out(1)=-RHS(1);  
  elseif (field.left.type=='V')
    out(1)=(M(1,1)+M(1,2))*field.internal(1)-RHS(1);  
  end
  
  if (field.right.type=='G' | field.right.type=='BP')
    out(end)=RHS(end);  
  elseif (field.right.type=='V')
    out(end)=(M(end,end)+M(end,end-1))*field.internal(end)-RHS(end);  
  end

end
