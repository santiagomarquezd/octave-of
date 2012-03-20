% Solves an 1D Mixture Model problem like OpenFOAM

%function mixtureSolver

page_screen_output(0);

  %tubopasante
  %tubopasanteTrivial
  %tubopasanteSaltoAlphag
  %tubopasanteSaltoV
  %tubopasanteMomCteVmayorInit
  %Garcia_Cascales_Phase_separation
  tubopasanteSaltoVmSaltoAlphag

% **************************** MAIN PROGRAM ***************************

% Gravity treatment terms
ghf = g*xF;
gh = g*xC;

% Set fields as 'old' states
rhom0=rhom;
alphag0=alphag;
U0=U;

% rhomPhi field initialization from U and rhom fields
rhomPhi=fvc_interpolate(U0, w, xC, xF).*fvc_interpolate(rhom0, w, xC, xF).*Sf;


% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

  % Mixture equation solution
  rhoEqn

  % Drift velocity calculation
  calcVdrp

  % UEqn
  UEqn
  
%if (step<2)  
  % PISO loop
  for corr=1:nCorr
    pEqnAlphaPISO
  end

  % Set fields as 'old' states
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
%end

end


