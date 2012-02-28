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


  % Volumetric flux from mass flux
  phiAlpha=rhomPhi./fvc_general_interpolate(rhom, xC, xF,-1,directionFluxU);

  %phiVdrp=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF).*Sf;
  directionFlux=fvc_interpolate(Vpq, w, xC, xF);
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiVdrp=fvc_general_interpolate(assign(Vpq,arrayToField(1-cp),'*'), xC, xF,1,directionFlux).*Sf;
  phiAlpha+=phiVdrp;

  % alphag0 values at interfaces with full upwind (direction given by phiAlpha)
  directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
  alphag0Int=fvc_general_interpolate(alphag0, xC, xF,-1,directionFluxAG0);

  %dt=-(1-alphag0.internal(end))*V(end)/...
      (phiAlpha(end).*alphag0Int(end)-phiAlpha(end-1).*alphag0Int(end-1));

  % Explicit solving for Alphag
  alphag.internal(1:end)=alphag0.internal(1:end)-dt./dx*(phiAlpha(2:end).*alphag0Int(2:end)-phiAlpha(1:end-1).*alphag0Int(1:end-1));

  % Artificial bounding
  % alphag=bound(alphag,'min',0);

  % BC setting
  alphag=setBC(alphag,rhom,xC,xF,g);
  %keyboard; pause;
  
  % Mixture density actualization
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);

  %keyboard; pause;

  % Set fields as 'old' states
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
  rhomPhi0=rhomPhi;
%end

end

end % Ends condition for temporal loop


