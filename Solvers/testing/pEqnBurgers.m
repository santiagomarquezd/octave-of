% pEqn.m
printf('PISO correction number:%d\n',corr)
%if (corr<4)
rUA.internal = 1./Aop(UEqnM, V);
rUA.left.type='V';
rUA.left.value=rUA.internal(1);
rUA.right.type='V';
rUA.right.value=rUA.internal(end);
rUA=setBC(rUA,constField(0,N),xC,xF,g);

rUAf=fvc_interpolate(rUA,w,xC,xF);

U.internal = rUA.internal.*Hop(UEqnM, UEqnRHS, U, V);
U=setBC(U,constField(0,N),xC,xF,g);

rhomPhi = 1/2*fvc_interpolate(rhom,w,xC,xF).*((fvc_interpolate(U,w,xC,xF).*Sf));
% fvc::ddtPhiCorr(rUA, rhom, U, rhomPhi)) correction NOT INCLUDED

rhomPhiU=rhomPhi;

% VERRRRRRRRRRRRRR ESTE TERMINO YA QUE CAMBIA DEBIDO A LA DIVISION DEL 
% RHS DE LA ECUACION DE MOMENTO POR RHOM
rhomPhi -= ghf.*fvc_snGrad(rhom,xC,xF).*rUAf.*Sf.*fvc_interpolate(rhom,w,xC,xF);

for nonOrth=1:nNonOrthCorr

    % Assembling (implicit terms)
    disp('Assembling p_rghEqn')
    [p_rghM, p_rghRHS]=fvm_laplacian(p_rgh,rUAf,xC,xF,Sf);

    % Calculating explicit terms (to RHS)
    ERHS=(fvc_ddt(rhom, rhom0, dt)+fvc_div_face(rhomPhi,V)).*V;

    %p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
    % **************** NOT IMPLEMENTED ********************

    % Solve
    disp('Solving for p_rgh')
    p_rgh.internal=p_rghM\(p_rghRHS+ERHS);  

end
if 1
rhomPhi -= flux(p_rghM,p_rghRHS,p_rgh);

% ****************************** PATCH *************************
%rhomPhi(1)=0;
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
p = assign(p_rgh,assign(rhom,arrayToField(gh),'*'),'+');

% **************** NOT IMPLEMENTED ********************


%  if (p_rgh.needReference())
%  {
%      p += dimensionedScalar
%      (
%          "p",
%          p.dimensions(),
%          pRefValue - getRefCellValue(p, pRefCell)
%      );
%      p_rgh = p - rhom*gh;
%  }

% ****************************************************

% Time marching for rhom Eqn. 23.4-1 from Fluent's User's Guide
% rhom0=rhom;
% It is a fixed point iteration over rho, so that rhom0 stays rhom
% from previous timestep

if 0
  rhoEqn
else
  alphaEqnBurgers
end

% Check continuity errors
% #include "mixtureContinuityErrs.H"       ---------------->>>>>>>> COMPLETAR

if 1
U.internal+=rUA.internal.*fvc_reconstruct((rhomPhi - rhomPhiU)./rUAf./fvc_interpolate(rhom,w,xC,xF),Sf);
U=setBC(U,constField(0,N),xC,xF,g);
end

end