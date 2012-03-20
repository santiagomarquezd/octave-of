% Solves an 1D Mixture Model problem like OpenFOAM

%function mixtureSolver

page_screen_output(0);

% Initialization
tubopasanteMomCteVmayorInit
%tubopasanteSaltoAlphagInit
%tubopasanteSaltoVmInit


% **************************** MAIN PROGRAM ***************************

% Gravity treatment terms
ghf = g*xF;
gh = g*xC;


% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

  % Mixture equation solution
  %rhoEqn

  
if 0
  % alphaEqn
  alphaEqn

  % UEqn
  UEqn
else
    
    
  if (strcmp(alphaVmScheme,'LxF'))
    [flux1,flux2]=no_Vr_cell_flux(U,alphag,rhom);
    [U,alphag]=LxF(U,alphag,flux1,flux2,dx,dt);
  elseif (strcmp(alphaVmScheme,'Rusanov'))
	% Local velocities at interfaces
	aeigens=aspeed(U,alphag,rhol,rhog,N);
    [flux1,flux2]=no_Vr_cell_flux(U,alphag,rhom);
    [U,alphag]=Rusanov(U,alphag,flux1,flux2,aeigens,dx,dt);
  elseif (strcmp(alphaVmScheme,'KT'))
    % Limiting
    % Sweby's fuction calculation
    phiAlphag=superbee(rvalue(alphag,1E-9)); 
    phiU=superbee(rvalue(U,1E-9));
    % Limited values calculation
    [ULimited]=limitedValues(U,phiU,dx,dt);
    [alphagLimited]=limitedValues(alphag,phiAlphag,dx,dt);
    % Local velocities at interfaces
    aeigens=aspeed(U,alphag,rhol,rhog,N);
    % Time advancement by Kurnanov & Tadmor's scheme
    [Udummy,alphag]=KT(U,alphag,ULimited,alphagLimited,@VmFlux,@alphaFlux,aeigens,dx,dt,0,rhol,rhog,1);
  end
  
  % Mixture density actualization
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
 
  % BC's actualization
  U=setBC(U,constField(0,N),xC,xF,g);
  alphag=setBC(alphag,rhom,xC,xF,g);
  %rhom=setBC(rhom,constField(0,N),xC,xF,g);
  
end

%end

end


