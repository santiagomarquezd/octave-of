% Upwind based in front velocities
% N umerical diffusivity
nu=mult*1/2*mean(abs(U.internal)+abs(Vpq.internal))*mean(dx)*ones(size(rhomPhi));

% nu overriden for test purposes
if 1
  nu=0.00705*ones(size(rhomPhi)); %00705
end

% Creates Alphag0 and Alphag as a copy of alpha.
Alphag0=alphag;
Alphag=Alphag0;

% Determination of Alpha auxiliary field
%Alphag0=assign(assign(alphag0,constField(rhog,N),'*'),rhom,'/');
Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
Alphag0=setBC(Alphag0,rhom,xC,xF,g);

% First part part of Alpha flux is given by the flux from PISO
phiAlpha=rhomPhi;

% Upwinding following wave velocities (sort of UADE stabilization)
if 0
  % Direction selected by inspection
  directionUp=(fvc_interpolate(alphag, w, xC, xF)>0.5)*(-2)+1;
  directionUp=directionUp(2:end-1,1);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,-1,directionUp).*Sf;
  phiAlpha+=phiVdrp;
    
  % Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
  Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionUp);
elseif 0
  % Direction selected by front velocity direction (discrete calculus)
  %phiCell=assign(assign(rhom,assign(U,assign(Vpq,arrayToField(1-cp),'*'),'+'),'*'),Alphag0,'*');
  phiCell=assign(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'),Alphag0,'*');
  % Front velocities calculation
  S=(phiCell.internal(2:end)-phiCell.internal(1:end-1)+1E-9)./(Alphag0.internal(2:end)-Alphag0.internal(1:end-1)+1E-9);
  directionUp=sign(S);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,-1,directionUp).*Sf;
  phiAlpha+=phiVdrp;
  % Alphag0 values at interfaces with full upwind (direction given by front velocity)
  Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionUp);
else
  % Direction selected by front velocity direction using eigenvalue of isolated alpha eqn
  alphaEqnIsolatedEigenvalues
  directionUp=fvc_interpolate(LL, w, xC, xF);
  directionUp=sign(directionUp);
  directionUp=directionUp(2:end-1,1);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,-1,directionUp).*Sf;
  phiAlpha+=phiVdrp;
  % Alphag0 values at interfaces with full upwind (direction given by front velocity)
  Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionUp);
end

  % Solve
if (fullVerbose==1)
  disp('Explicit solving of Alphag')
end

if 1
  % Forces boundary fluxes to be zero
  phiAlpha(1)=0;
  phiAlpha(end)=0;
end	

% Explicit solving for Alphag
Alphag.internal(1:end)=Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx*(phiAlpha(2:end).*...
  			 Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1))./rhom.internal(1:end);

Alphag=setBC(Alphag,rhom,xC,xF,g);

% Mixture density actualization
rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);
% BC overriden (maybe rhom has to have ZG BC at top and bottom)
rhom.left.value=rhom.internal(1);
rhom.right.value=rhom.internal(end);
rhom=setBC(rhom,constField(0,N),xC,xF,g);

% alpha actualization
alphag.internal=rhom.internal.*Alphag.internal./rhog;
alphag=setBC(alphag,rhom,xC,xF,g);

%load densities.dat
