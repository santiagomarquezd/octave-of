% Solves an 1D Mixture Model problem in diagonalized form

%function mixtureSolver

page_screen_output(0);

% Initialization
tubopasanteMomCteVmayorInit
%tubopasanteSaltoAlphagInit
%tubopasanteSaltoVmInit


% **************************** MAIN PROGRAM ***************************


% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)
  
  % Transformation from primitive variables to characterist ones
  [char1,char2]=primToChar(U,alphag,@eigens,rhol,rhog)
  
  % Once variables have been transformed advection is done.

end


