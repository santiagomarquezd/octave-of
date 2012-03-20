% UEqn.m
% Solve the Momentum equation
% Eqn. 23.4-4 from Fluent's User's Guide

% *************************************
% No viscous term is taken into account
% *************************************


    
% Assembling (implicit terms)
disp('Assembling UEqn')
%[M, RHS] = fvm_ddt(rhom,U0,V,dt,1);
% Lagged rho version (???????????????)
[M, RHS] = fvm_ddt(constField(1,N),constField(1,N),U0,V,dt,1);
[A, ARHS]= fvm_div_flux_cell(rhomPhi./fvc_interpolate(rhom,w,xC,xF),U0,xC,xF,w,1);

% Calculating explicit terms (to RHS)
% cp updating
cp=alphag.internal.*rhog./rhom.internal;
arg=assign(assign(rhom,arrayToField((1-cp).*cp),'*'),assign(Vpq,constField(2,N),'^'),'*');
ERHS=-fvc_div_cell(arg, w, xC, xF, Sf, V).*V;
% Pressure terms are calculated apart in order to correctly calculate the H operator
FRHS=fvc_reconstruct((-ghf.*fvc_snGrad(rhom,xC,xF)-fvc_snGrad(p_rgh,xC,xF))./fvc_interpolate(rhom,w,xC,xF).*Sf,Sf).*V;


% Final assembling
UEqnM=M+A;
UEqnRHS=RHS+ARHS+ERHS;

% Solve
disp('Solving for U')
U.internal=UEqnM\(UEqnRHS+FRHS); 
