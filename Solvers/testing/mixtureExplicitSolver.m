% Solves an 1D Mixture Model using the alpha equation
% and algebraic momentum equation

clear all;
page_screen_output(0);

  Garcia_Cascales_Phase_separation

% **************************** MAIN PROGRAM ***************************

% Set fields as 'old' states
rhom0=rhom;
alphag0=alphag;
U0=U;

% rhomPhi field initialization from U and rhom fields
% Creation of flux direction
directionFlux=fvc_interpolate(U0, w, xC, xF);
directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
rhomPhi=fvc_interpolate(U0, w, xC, xF).*fvc_general_interpolate(rhom0, xC, xF,-1,directionFlux).*Sf;

% Set fields as 'old' states
rhomPhi0=rhomPhi;

if 1 % Enables temporal loop

% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

  % Drift velocity calculation
  calcVdrp
  % calcVdrpSimple

  % UEqn
  UEqnSimple

  % rhomPhi field calculation from U and rhom fields
  % Creation of flux direction
  directionFluxU=fvc_interpolate(U, w, xC, xF);
  directionFluxU=sign(directionFluxU(2:end-1,1)+1E-9);
  %rhomPhi=fvc_interpolate(U, w, xC, xF).*fvc_general_interpolate(rhom, xC, xF,-1,directionFluxU).*Sf;
  rhomPhi=fvc_general_interpolate(U, xC, xF,-1,directionFluxU).*fvc_general_interpolate(rhom, xC, xF,-1,directionFluxU).*Sf;

  % alphaEqn explicit calculation
  % Numerical diffusivity
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

  phiAlpha=rhomPhi;

  %phiVdrp=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF).*Sf;
  directionFlux=fvc_interpolate(Vpq, w, xC, xF);
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;
  phiAlpha+=phiVdrp;

  % Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
  directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
  Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionFluxAG0);

  %dt=-(1*rhog/1-Alphag0.internal(end))*V(end)/...
  %    (phiAlpha(end).*Alphag0Int(end)-phiAlpha(end-1).*Alphag0Int(end-1));

  % Explicit solving for Alphag
  Alphag.internal(1:end)=Alphag0.internal(1:end)-dt./dx*(phiAlpha(2:end).*Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1));

  % Artificial bounding
  %Alphag=bound(Alphag,'min',0);

  % BC setting
  Alphag=setBC(Alphag,rhom,xC,xF,g);
  %keyboard; pause;

if 1 
  % Mixture density actualization
  rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);  %ORIGINAL
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);

  % alpha actualization
  alphag.internal=rhom.internal.*Alphag.internal./rhog;   % ORIGINAL
  alphag=setBC(alphag,rhom,xC,xF,g);

else

  % alpha actualization
  alphag.internal=rhom0.internal.*Alphag.internal./rhog;
  alphag=setBC(alphag,rhom,xC,xF,g);

  % Mixture density actualization
  rhom.internal=rhog.*alphag.internal./Alphag.internal;
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);

end


  %keyboard; pause;

  % Set fields as 'old' states
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
  rhomPhi0=rhomPhi;
%end

end

end % Ends condition for temporal loop


