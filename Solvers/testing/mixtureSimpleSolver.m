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

  % Mixture equation solution
  % rhoEqn
 
  % Drift velocity calculation
  calcVdrp
  % calcVdrpSimple

  % UEqn
  UEqnSimple

  % alphaEqn
  %alphaEqn
  alphaEqnIshiiSimple

  % Set fields as 'old' states
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
  rhomPhi0=rhomPhi;
%end

end

end % Ends condition for temporal loop


