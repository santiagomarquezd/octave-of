% UEqnSimple.m
% Solve the Momentum equation
% via alebraic formulas

U.internal=alphag.internal.*(rhog./rhom.internal-1).*Vpq.internal; 
U=setBC(U,rhom,xC,xF,g);



