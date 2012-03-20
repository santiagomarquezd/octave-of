% Solves an 1D Mixture Model problem like OpenFOAM

%function mixtureSolver

page_screen_output(0);

readinitfile=1

  %tubopasante
  %tubopasanteTrivial
  %tubopasanteSaltoAlphag
  %tubopasanteSaltoV
  tubopasanteMomCteVmayorInit

% **************************** MAIN PROGRAM ***************************

% Gravity treatment terms
ghf = g*xF;
gh = g*xC;

% Set fields as 'old' states
rhom0=rhom;
alphag0=alphag;
U0=U;

% phi field initialization from U
%phi=1/2*fvc_interpolate(U0, w, xC, xF).*Sf;

% rhomPhi field initialization from U and rhom fields
rhomPhi=1/2*fvc_interpolate(U0, w, xC, xF).*fvc_interpolate(rhom0, w, xC, xF).*Sf;



% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

  % Mixture equation solution
  % rhoEqn



  % Drift velocity calculation
  calcVdrp
  
  % UEqn
  UEqnBurgers
 
  % alphaEqn
  alphaEqnBurgers
 
  %PISO loop
    for corr=1:nCorr
      pEqnBurgers
    end

%if (step<0)
  % Set fields as 'old' states
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
%end
  
end


