% Solves an 1D Mixture Model problem like OpenFOAM

%function mixtureSolver

clear all;
page_screen_output(0);

  %tubopasante
  %tubopasanteTrivial
  %tubopasanteSaltoAlphag
  %tubopasanteSaltoV
  %tubopasanteMomCteVmayorInit
  Garcia_Cascales_Phase_separation

% **************************** MAIN PROGRAM ***************************

% Gravity treatment terms
ghf = g*xF;
gh = g*xC;

%  pp=(-(1-xC).*rhom.internal*g)-rhom.internal.*gh;
%  p_rgh.internal=pp;
%  p_rgh.right.type='V';
%  p_rgh.right.value=p_rgh.internal(end);
%  p_rgh=setBC(p_rgh,rhom,xC,xF,g);


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

if (initia!=1)
  load("dump.dat");
end

if 1 % Enables temporal loop

% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

  % Mixture equation solution
  rhoEqn
 
  % Drift velocity calculation
  calcVdrp

  % UEqn
  UEqn

  % alphaEqn
  %alphaEqn
  alphaEqnIshii

  %PISO loop
  if 1  
    for corr=1:nCorr
      pEqn
    end
  end

%if (step<timesteps)

  % Set fields as 'old' states
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
  rhomPhi0=rhomPhi;
%end

end

end % Ends condition for temporal loop

save -binary dump.dat rhom0 rhom alphag0 alphag U0 U rhomPhi0 rhomPhi Alphag0 Alphag

%end % End function mixtureSolver
