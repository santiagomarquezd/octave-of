% pEqn.m
if (fullVerbose==1)
  printf('PISO correction number:%d\n',corr)
end

%if (corr<4)
rUA.internal = 1./Aop(UEqnM, V);
rUA.left.type='V';
rUA.left.value=rUA.internal(1);
rUA.right.type='V';
rUA.right.value=rUA.internal(end);
rUA=setBC(rUA,constField(0,N),xC,xF,g);

rUAf=fvc_interpolate(rhom,w,xC,xF).*fvc_interpolate(rUA,w,xC,xF);


U.internal = rUA.internal.*Hop(UEqnM, UEqnRHS, U, V);
U=setBC(U,constField(0,N),xC,xF,g);

phiC = fvc_ddtPhiCorrection(U0, rhomPhi0, rUA, rUAf, rhom0, w, xC, xF, Sf, dt);
%keyboard; pause;

rhomPhi = fvc_interpolate(rhom,w,xC,xF).*((fvc_interpolate(U,w,xC,xF).*Sf)+phiC);

%if (step<timesteps || corr<stopCorr)

rhomPhiU=rhomPhi;

rhomPhi -= ghf.*fvc_snGrad(rhom,xC,xF).*rUAf.*Sf;

% This line has been added but it isn't in original code
rhomPhi(end)=0;

%if 0

for nonOrth=1:nNonOrthCorr

    % Assembling (implicit terms)
    if (fullVerbose==1)
      disp('Assembling p_rghEqn')
    end
    [p_rghM, p_rghRHS]=fvm_laplacian(p_rgh,rUAf,xC,xF,Sf);

    % Calculating explicit terms (to RHS)
    ERHS=(fvc_ddt(rhom, rhom0, dt)+fvc_div_face(rhomPhi,V)).*V;

    %p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
    % **************** NOT IMPLEMENTED ********************

    % Solve
    if (fullVerbose==1)
      disp('Solving for p_rgh')
    end
    p_rgh.internal=p_rghM\(p_rghRHS+ERHS);  

%      if (nonOrth == nNonOrthCorr)
%      {
%          rhomPhi -= p_rghEqn.flux();
%      }
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

rhoEqn

% Check continuity errors
% #include "mixtureContinuityErrs.H"       ---------------->>>>>>>> COMPLETAR

if 1
U.internal+=rUA.internal.*fvc_reconstruct((rhomPhi - rhomPhiU)./rUAf,Sf);
U=setBC(U,constField(0,N),xC,xF,g);
end

end % End of conditional part of code

%end % Conditional for timestep and correction