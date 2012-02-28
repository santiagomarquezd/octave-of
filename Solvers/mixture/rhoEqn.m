%rhoEqn.m

% Assembling
if (fullVerbose==1)
  disp('Assembling rhomEqn')
end
[rhoMat, RHS] = fvm_ddt(constField(1,N),constField(1,N),rhom0,V,dt,1);
%[A, ARHS]= fvm_div_flux_cell(rhomPhi,rhom0,xC,xF,w,1);
ARHS=-fvc_div_face(rhomPhi,V).*V ; 
%  M=M+A;
%  RHS=RHS+ARHS;  
 
% Solve
if (fullVerbose==1)
  disp('Solving for rhom')
end
rhom.internal=rhoMat\(RHS+ARHS);

% BC overriden (maybe rhom has to have ZG BC at top and bottom)
%  rhom.left.value=rhom.internal(1);
%  rhom.right.value=rhom.internal(end);
%  rhom=setBC(rhom,constField(0,N),xC,xF,g);