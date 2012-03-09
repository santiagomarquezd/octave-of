function [out]=setBC(in,rho,xC,xF,g)
  % Evaluate BC's 
  %
  % [out]=setBC(in,rho,xC,xF,g)
  %
  % out: field with evaluated BC's
  % in: field without evaluated BC's
  % rho: problem density field for Buoyant Pressure BC's
  % xC: cell centres
  % xF: face centres
  % g: problem gravity for Buoyant Pressure BC's
  
  % General copy
  out=in;
  
  % Left Boundary conditions
  if (in.left.type=='V')
	out.left.setvalue=in.left.value;
  elseif (in.left.type=='G')
	out.left.setvalue=in.internal(1)-in.left.gradient*(xC(1)-xF(1));
  elseif (in.left.type=='BP')
	rhoSnGrad=fvc_snGrad(rho,xC,xF);
	out.left.gradient = -rhoSnGrad(1)*g*xF(1);
	out.left.setvalue=in.internal(1) + out.left.gradient*abs(xC(1)-xF(1));
  end
  
  % Right Boundary conditions
  if (in.right.type=='V')
	out.right.setvalue=in.right.value;
  elseif (in.right.type=='G')
	out.right.setvalue=in.internal(end)+in.right.gradient*(xF(end)-xC(end));
  elseif (in.right.type=='BP')
	rhoSnGrad=fvc_snGrad(rho,xC,xF);
	out.right.gradient = -rhoSnGrad(end)*g*xF(end);
	out.right.setvalue=in.internal(end) + out.right.gradient*abs(-xF(end)-xC(end));
  end	

end