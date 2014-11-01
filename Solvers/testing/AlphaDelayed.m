% ----------------------------------------------------------------------
% Explicit solution of Alpha by UADE method with Delayed rhom (Gastaldo)

% alphaEqn
phiAlpha=rhomPhi;

% Final flux, PISO part
fluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF);

% UADE stabilization with full upwind/downwind
directionFlux=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF)+phiAlpha;
directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;

phiAlpha+=phiVdrp;

% Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionFluxAG0);

% Final flux, relative velocity part
fluxAlpha=fluxAlpha+phiVdrp.*Alphag0Int;

% Solve
if (fullVerbose==1)
  disp('Explicit solving of Alphag')
end	

if 0
      % Impermeable walls test
      phiAlpha(1)=0;
      phiAlpha(end)=0;
end

% Explicit solving for AlphagRhom
if 1
  % Overwrites the flux with a new one where Alpha is ever Downwind/Upwind
  fluxAlpha=phiAlpha.*Alphag0Int;
end

Alphag.internal(1:end)=Alphag0.internal(1:end).*rhom0bak.internal(1:end)./rhom.internal(1:end)-dt./dx*(fluxAlpha(2:end)-...
			fluxAlpha(1:end-1))./rhom.internal(1:end);

Alphag=setBC(Alphag,rhom,xC,xF,g);

% -----------------------------------------
% Explicit solution of Alpha by KTcFlux method

% New version based in K&T scheme, High Resolution schemes can be activated (testing)
% ****************************************************************************  
% IMPORTANT: Flux calculation and stabilization is based directly on alphag
% ****************************************************************************

% Limiting
% Sweby's fuction calculation
% phiAlphag=superbee(rvalue(alphag0,1E-9));
% phiAlphag=vanLeer(rvalue(alphag0,1E-9));
 phiAlphag=vanLeer(rvalue(alphag0,1E-9))*0; % Constant values by cells (mimiking Rusanov?)

% Limited values calculation
[alphagLimited]=limitedValues(alphag0,phiAlphag,dx,dt);

% Here it'd be used eigenvalues calculated as dF/dA like in aspeesIsolatedAlphaEqnPhiVmRhoT.m
% it has into account the frozen values of Vm*rhom given in rhomPhi and the eventual frozen value of 
% rhom prediction rhomT. Both of these values, rhomPhi and rhonT (rhom) are given.  
 [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(alphag0,rhol,rhog,V0,aexp);    
% [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVmGastaldo(alphag,rhol,rhog,V0); % Gastaldo's example
% Stabilization flux has to be zero at boundaries
a_j_minus_half(1)=0;
a_j_plus_half(end)=0;
% Flatten for KT function
aeigens=[a_j_minus_half; a_j_plus_half(end)];

% Alpha at faces
AlphagInt=fvc_interpolate(Alphag0, w, xC, xF);

% rhom0 at faces
rhom0bakInt=fvc_interpolate(rhom0bak, w, xC, xF);

% Time advancement by Kurnanov & Tadmor's scheme
[AlphagRF]=KTcFlux(Alphag0,alphagLimited,@AEqnNoVmFluxFlat,aeigens,rhomPhi,AlphagInt,dx,dt,rhom,rhom0bak,rhom0bakInt,V0,rhol,rhog,aexp);


% -----------------------------------------
% Seleccion of method of temporal actulization

disp('A antes de actualizar')
sum(Alphag0.internal)
disp('A despues de actualizar (UADE)')
sum(Alphag.internal)
disp('A despues de actualizar (Rusanov)')
sum(AlphagRF.internal)

% Uses UADE solution (0) or Rusanov (1)
if 1
  Alphag=AlphagRF;
end

% -----------------------------------------
% Temporal actualization of remaining variables (rhom and alphag)

% Mixture density actualization
rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);
% BC overriden (maybe rhom has to have ZG BC at top and bottom)
rhom.left.value=rhom.internal(1);
rhom.right.value=rhom.internal(end);
rhom=setBC(rhom,constField(0,N),xC,xF,g);

% alpha actualization
alphag.internal=rhom.internal.*Alphag.internal./rhog;
alphag=setBC(alphag,rhom,xC,xF,g);

rhomSave=rhom;


