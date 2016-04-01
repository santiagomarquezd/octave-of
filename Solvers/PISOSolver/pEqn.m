% pEqn.m
if (fullVerbose == 1)
  printf('PISO correction number:%d\n', corr)
end

rUA.internal = 1./Aop(UEqnM, V);
rUA.left.type = 'V';
rUA.left.value = rUA.internal(1);
rUA.right.type = 'V';
rUA.right.value = rUA.internal(end);
rUA = setBC(rUA, constField(0,N), xC, xF, g);

rUAf = fvc_interpolate(rUA, w, xC, xF);

U.internal = rUA.internal.*Hop(UEqnM, UEqnRHS, U, V);
U = setBC(U, constField(0,N), xC, xF, g);

% Set to zero
phiC = fvc_ddtPhiCorrection(U0, phi0, rUA, rUAf, ... 
rho0, w, xC, xF, Sf, dt)*0;

phi = fvc_interpolate(rho, w, xC, xF).* ...
 ((fvc_interpolate(U, w, xC, xF).*Sf) + phiC);

phi += g.*fvc_interpolate(rho, w, xC, xF).*Sf.*rUAf;

% Assembling (implicit terms)
if (fullVerbose == 1)
    disp('Assembling p')
end

% p needs to be fixed
[pM, pRHS] = fvm_laplacian(p, rUAf, xC, xF, Sf);

% Calculating explicit terms (to RHS)
ERHS = fvc_div_face(phi, V).*V;

% Solve
if (fullVerbose == 1)
    disp('Solving for p')
end

p.internal = pM\(pRHS + ERHS);  

% Flux correction
phi -= flux(pM, pRHS, p);

U.internal += rUA.internal.*(g - fvc_grad(p, w, xC, xF, Sf, V));
U = setBC(U, constField(0,N), xC, xF, g);

