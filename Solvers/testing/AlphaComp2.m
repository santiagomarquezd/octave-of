% Artificial flux
if 1
  Utmp=assign(Vpq,assign(alphag,assign(assign(constField(rhog,N),rhom,'/'),constField(1,N),'-'),'*'),'*');
  if 1
    % Free upper boundary velocity
    Utmp.right.type='G';
    Utmp.right.gradient=0;
  end
  % BC's setting
  Utmp=setBC(Utmp,constField(0,N),xC,xF,0);
  U=Utmp; % Overwrite the momentum predictor
  % Flux updating
  %rhomPhi=fvc_interpolate(assign(rhom,U,'*'), w, xC, xF)
  rhomPhi=fvc_interpolate(rhom, w, xC, xF).*fvc_interpolate(Utmp, w, xC, xF);
end

% -----------------------------------------
% Explicit solution of Alpha by UADE method

% To avoid oscillations
%if (rhomPhi(end)<0) rhomPhi=0; end

% alphaEqn
phiAlpha=rhomPhi;

% Final flux, PISO part
fluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF).*1;

% UADE stabilization with full upwind/downwind
directionFlux=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF)+phiAlpha;
directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;

phiAlpha+=phiVdrp;

if 0
  % Impermeable walls test
  phiAlpha(1)=0;
  phiAlpha(end)=0;
end
  
% Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionFluxAG0);

% Final flux, relative velocity part
fluxAlpha=fluxAlpha+phiVdrp.*Alphag0Int;

% Solve
if (fullVerbose==1)
  disp('Explicit solving of AlphagRhom')
end

if 1
  % Overwrites the flux with a new one where Alpha is ever Downwind/Upwind
  fluxAlpha=phiAlpha.*Alphag0Int;
end
	

% Explicit solving for AlphagRhom
phiAlpha=phiAlpha.*Alphag0Int;
AlphagRhom.internal(1:end)=AlphagRhom0.internal(1:end)-dt./dx*(phiAlpha(2:end)-phiAlpha(1:end-1));
AlphagRhom=setBC(AlphagRhom,rhom,xC,xF,g);

% -----------------------------------------
% Explicit solution of Alpha by Rusanov method

% Face maximum for spectral radius specially tailored for this problem
if 0
  % U-alphag block
  [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag0,rhol,rhog,V0);
else
  % Isolated alpha equation with explicit relation for Vm (not calculating dF/dalphag)
  [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(alphag,rhol,rhog,V0);
end

if 1
  % Stabilization flux has to be zero at boundaries too (impermeable wall)
  a_j_minus_half(1)=0;
  %a_j_plus_half(end)=0;  
end

% Quasi-conservative, rhomPhi is used instead of creating it
fluxVdrpRF=assign(AlphagRhom0,assign(Vpq,arrayToField(1-cp),'*'),'*');
fluxVdrpRF=setBC(fluxVdrpRF,constField(0,N),xC,xF,0);

if 0
  %Ensure no fluxes at boundaries
  fluxVdrpRF.left.type='V';
  fluxVdrpRF.left.value=0;
  fluxVdrpRF.right.type='V';
  fluxVdrpRF.right.value=0;
  fluxVdrpRF=setBC(fluxAlpha,constField(0,N),xC,xF,0);
end

% Manual solution taking values at faces
% Flux calculation
% Internal faces
fluxVdrpRFf=zeros(N+1,1);
% Standard central difference  
fluxVdrpRFf(2:end-1)=(fluxVdrpRF.internal(1:end-1)+fluxVdrpRF.internal(2:end))/2;

% Stabilization
% Final flux, internal faces
% Here it'd be used eigenvalues calculated as dF/dA like in aspeesIsolatedAlphaEqnPhiVmRhoT.m
% it has into account the frozen values of Vm*rhom given in rhomPhi and the eventual frozen value of 
% rhom prediction rhomT. Both of these values, rhomPhi and rhonT (rhom) are given.  
% *********************************************
% The eigenvalues have to be multiplied by rhog
% due the flux has units of V*rhog
% ********************************************* 
fluxVdrpRFf(2:end-1)=fluxVdrpRFf(2:end-1)-rhog.*a_j_plus_half(1:end-1).*(Alphag0.internal(2:end)-Alphag0.internal(1:end-1))/2;
% Boundary faces
fluxVdrpRFf(1)=fluxVdrpRF.left.setvalue+rhog*a_j_minus_half(1).*(Alphag0.internal(1)-Alphag0.left.setvalue);
fluxVdrpRFf(end)=fluxVdrpRF.right.setvalue+rhog*a_j_plus_half(1).*(Alphag0.right.setvalue-Alphag0.internal(end));
% Flux from PISO is added multiplied by the corresponding Alphag0 values
Alphag0f=fvc_interpolate(Alphag0, w, xC, xF);

if 1
  % Use central difference for Alphag0f
  fluxAlphaRFf=fluxVdrpRFf+rhomPhi.*Alphag0f;
elseif 1
  % Use upwinding for Alphag0f, this value is given by Alphag0Int calculated in UADE's section
  fluxAlphaRFf=fluxVdrpRFf+rhomPhi.*Alphag0Int;
elseif 0
  % Flux for mixture velocity is not taken from PISO
  fluxAlphaRFf=fluxVdrpRFf+fvc_interpolate(assign(assign(rhom,U,'*'),Alphag0,'*'), w, xC, xF);
else
  % Use semi-analytic flux
  fluxAlphaRFf=fvc_interpolate(fluxAlphaRFAnalytic, w, xC, xF);
end

% Time advancement
% Memory allocation
AlphagRhomRF=AlphagRhom0;
% Calculation
AlphagRhomRF.internal(1:end)=AlphagRhom0.internal-dt./dx.*(fluxAlphaRFf(2:end)-fluxAlphaRFf(1:end-1));
AlphagRhomRF=setBC(AlphagRhomRF,rhom,xC,xF,g);


% -----------------------------------------
% Temporal actualization of remaining variables (rhom and alphag)

disp('ARhom (alphag*rhog) antes de actualizar')
sum(AlphagRhom0.internal)
disp('ARhom (alphag*rhog) despues de actualizar (UADE)')
sum(AlphagRhom.internal)
disp('ARhom (alphag*rhog) despues de actualizar (Rusanov)')
sum(AlphagRhomRF.internal)

% Uses UADE solution (0) or Rusanov (1)
if 0
  AlphagRhom=AlphagRhomRF;
end

% alpha actualization
alphag.internal=AlphagRhom.internal./rhog;
alphag=setBC(alphag,rhom,xC,xF,g);

% Mixture density actualization
rhom.internal=alphag.internal*rhog+(1-alphag.internal).*rhol;
% BC overriden (maybe rhom has to have ZG BC at top and bottom)
rhom.left.value=rhom.internal(1);
rhom.right.value=rhom.internal(end);
rhom=setBC(rhom,constField(0,N),xC,xF,g);

% Determination of actualized Alpha auxiliary field (consistent with alphaEqn solving
% not including rhom PISO changes)
Alphag.internal=alphag.internal.*rhog./rhom.internal;
Alphag=setBC(Alphag,rhom0,xC,xF,g);

