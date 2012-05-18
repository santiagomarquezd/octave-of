% Gastaldo's alphaEqn for post PISO correction of alpha and rhom

% Numerical diffusivity
nu=mult*1/2*mean(abs(U.internal)+abs(Vpq.internal))*mean(dx)*ones(size(rhomPhi));

% alphaEqn
if 1
  % UADE stabilization
  directionFlux=fvc_interpolate(Vpq, w, xC, xF);
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiAlpha=fvc_general_interpolate(assign(Vpq,arrayToField(1-cp),'*'), xC, xF,1,directionFlux).*Sf;
else
  % FOAM original
  phiAlpha=fvc_interpolate(assign(Vpq,arrayToField(1-cp),'*'), w, xC, xF).*Sf;
end

% Assembling
if (fullVerbose==1)
  disp('Assembling alphaGEqn')
end
[ddtMat, ddtRHS] = fvm_ddt(constField(1,N),constField(1,N),alphag0,V,dt,1);
% alphag has been previously actualized by Vm in in Pre stage
[divMat, divRHS]= fvm_div_flux_cell(phiAlpha,alphag,xC,xF,w,alphaDiv);
[lapMat, lapRHS]=fvm_laplacian(alphag0,nu,xC,xF,Sf);
  
%  M=M+A+B;
%  RHS=RHS+ARHS+BRHS;
M=ddtMat+divMat+lapMat;
RHS=ddtRHS+divRHS+lapRHS;
 
% Solve
if (fullVerbose==1)
  disp('Solving for alphag (post-advancing)')
end
alphag.internal=M\RHS;

% Secondary phase bounding
if 0
  alphag=bound(alphag,'min',0);
  alphag=bound(alphag,'max',1);
  alphag=setBC(alphag,rhom,xC,xF,g);
end

if 1
  % Mixture density actualization
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);
end
