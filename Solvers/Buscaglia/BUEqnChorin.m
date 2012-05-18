% UEqn.m
% Solve the Momentum equation
% Eqn. 23.4-4 from Fluent's User's Guide

% *************************************
% No viscous term is taken into account
% *************************************
    
% Assembling (implicit terms)
if (fullVerbose==1)
  disp('Assembling UEqn')
end
%[M, RHS] = fvm_ddt(rhom,U0,V,dt,1);
% Lagged rho version (???????????????)
[ddtM, ddtRHS] = fvm_ddt(rhom,rhom0,U0,V,dt,1);
[convM, convRHS]= fvm_div_flux_cell(rhomPhi,U0,xC,xF,w,1);

% Calculating explicit terms (to RHS)
% cp updating
cp=alphag.internal.*rhog./rhom.internal;
arg=assign(assign(rhom,arrayToField((1-cp).*cp),'*'),assign(Vpq,constField(2,N),'^'),'*');
driftRHS=-fvc_div_cell(arg, w, xC, xF, Sf, V).*V;
% No pressure term but gravity
volRHS=g.*rhom.internal.*V;

% Final assembling
UEqnM=ddtM+convM;
UEqnRHS=ddtRHS+convRHS+driftRHS;

% Solve
if (fullVerbose==1)
  disp('Solving for U')
end
U.internal=UEqnM\(UEqnRHS+volRHS);
U=setBC(U,constField(0,N),xC,xF,g);

% U from momentum predictor is stored for debugging purposes
UmomPred=U;

%keyboard; pause;

% U bounding
if 0
U=bound(U,'min',0);
U=bound(U,'max',1);
U=setBC(U,rhom,xC,xF,g);
end

%keyboard; pause;
%hold on;plot(xC,U.internal,'g*-');