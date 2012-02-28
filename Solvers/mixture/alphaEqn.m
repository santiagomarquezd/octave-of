% Numerical diffusivity
nu=mult*1/2*mean(abs(U.internal)+abs(Vpq.internal))*mean(dx)*ones(size(rhomPhi));

% alphaEqn
% Implements Eqn. 23.4-17 from Fluent's User's Guide
phiAlpha=rhomPhi./fvc_interpolate(rhom, w, xC, xF);
phiVdrp=fvc_interpolate(assign(Vpq,arrayToField(1-cp),'*'), w, xC, xF).*Sf;
phiAlpha+=phiVdrp;

% Assembling
disp('Assembling alphaGEqn')
[ddtMat, ddtRHS] = fvm_ddt(constField(1,N),constField(1,N),alphag0,V,dt,1);
[divMat, divRHS]= fvm_div_flux_cell(phiAlpha,alphag0,xC,xF,w,alphaDiv);
[lapMat, lapRHS]=fvm_laplacian(alphag0,nu,xC,xF,Sf);
  
%  M=M+A+B;
%  RHS=RHS+ARHS+BRHS;
M=ddtMat+divMat+lapMat;
RHS=ddtRHS+divRHS+lapRHS;
 
% Solve
disp('Solving for alphag')
alphag.internal=M\RHS;

% Secondary phase bounding
if 1
alphag=bound(alphag,'min',0);
alphag=bound(alphag,'max',1);
alphag=setBC(alphag,rhom,xC,xF,g);
end
%if i<2

% Mixture density actualization
rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
% BC overriden (maybe rhom has to have ZG BC at top and bottom)
rhom.left.value=rhom.internal(1);
rhom.right.value=rhom.internal(end);
rhom=setBC(rhom,constField(0,N),xC,xF,g);

%end