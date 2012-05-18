% UEqnVisco.m
% Solve the Momentum equation
% Eqn. 23.4-4 from Fluent's User's Guide

% Mixture dynamic viscosity
mu=nug*rhog*alphag.internal+(1-alphag.internal)*nul*rhol;

% Numerical diffusivity
muEst=multU*1/2*mean(abs(U.internal))*mean(dx)*ones(size(U)).*rhom.internal;

% Final dynamic viscosity
mu=mu+muEst;

% Interpolate in faces
mu=arrayToField(mu);
muf=fvc_interpolate(mu, w, xC, xF);
   
% Assembling (implicit terms)
if (fullVerbose==1)
  disp('Assembling UEqn')
end
%[M, RHS] = fvm_ddt(rhom,U0,V,dt,1);
% Lagged rho version (???????????????)
[ddtM, ddtRHS] = fvm_ddt(rhom,rhom0,U0,V,dt,1);
[convM, convRHS]= fvm_div_flux_cell(rhomPhi,U0,xC,xF,w,1);
[lapUM, lapURHS]= fvm_laplacian(U0,muf,xC,xF,Sf);

% Calculating explicit terms (to RHS)
% cp updating
cp=alphag.internal.*rhog./rhom.internal;
arg=assign(assign(rhom,arrayToField((1-cp).*cp),'*'),assign(Vpq,constField(2,N),'^'),'*');
driftRHS=-fvc_div_cell(arg, w, xC, xF, Sf, V).*V;
% - (fvc::grad(U) & fvc::grad(turbulence->muEff()))
viscoRHS=fvc_grad(U0, w, xC, xF, Sf, V).*fvc_grad(mu, w, xC, xF, Sf, V);
% Pressure terms are calculated apart in order to correctly calculate the H operator
volRHS=fvc_reconstruct((-ghf.*fvc_snGrad(rhom,xC,xF)-fvc_snGrad(p_rgh,xC,xF)).*Sf,Sf).*V;

% Final assembling
UEqnM=ddtM+convM-lapUM;
UEqnRHS=ddtRHS+convRHS-lapURHS+driftRHS; %+viscoRHS

% Solve
if (fullVerbose==1)
  disp('Solving for U')
end
U.internal=UEqnM\(UEqnRHS+volRHS); 

% U from momentum predictor is stored for debugging purposes
UmomPred=U;

% U bounding
if 0
U=bound(U,'min',0);
U=bound(U,'max',1);
U=setBC(U,rhom,xC,xF,g);
end

%keyboard; pause;