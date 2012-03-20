% alphaEqn
% Implements Eqn. 23.4-17 from Fluent's User's Guide

phiVdrp=fvc_interpolate(assign(Vpq,arrayToField(1-cp),'*'), w, xC, xF).*Sf;
phiAlpha=phiVdrp+rhomPhi./fvc_interpolate(rhom,w,xC,xF);

% Assembling
disp('Assembling alphaGEqn')
[M, RHS] = fvm_ddt(constField(1,N),constField(1,N),alphag0,V,dt,1);
[A, ARHS]= fvm_div_flux_cell(phiAlpha,alphag0,xC,xF,w,1);
%[B, BRHS]= fvm_div_flux_cell(phiVdrp,alphag0,xC,xF,w,1);
  
%  M=M+A+B;
%  RHS=RHS+ARHS+BRHS;
M=M+A;
RHS=RHS+ARHS;
 
% Solve
disp('Solving for alphag')
alphag.internal=M\RHS;

% Secondary phase bounding
if 0
alphag=bound(alphag,'min',0);
alphag=bound(alphag,'max',1);
alphag=setBC(alphag,rhom,xC,xF,g);
end
%if i<2

% Mixture density actualization
rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
rhom=setBC(rhom,constField(0,N),xC,xF,g);

%end