% UEqn.m
% Solve the Momentum equation


% *************************************
% No viscous term is taken into account
% *************************************
    
% Assembling (implicit terms)
if (fullVerbose==1)
    disp('Assembling UEqn')
end

[convM, convRHS] = fvm_div_flux_cell(phi, U0, xC, xF, w, 1);
[diffM, diffRHS] = fvm_laplacian(U, nuf, xC, xF, Sf);

% Calculating explicit terms (to RHS)

% Pressure terms are calculated apart in order
% to correctly calculate the H operator
volRHS = (g - fvc_grad(p, w, xC, xF, Sf, V)).*V;

% Final assembling
UEqnM = convM - diffM;
UEqnRHS = convRHS - diffRHS;

[UEqnM, UEqnRHS] = relaxM(UEqnM, UEqnRHS, U0, omegaU);

% Solve
if (fullVerbose == 1)
    disp('Solving for U')
end

U.internal = UEqnM\(UEqnRHS + volRHS); 

% U from momentum predictor is stored for debugging purposes
UmomPred = U;

%keyboard; pause;
%hold on;plot(xC,U.internal,'g*-');