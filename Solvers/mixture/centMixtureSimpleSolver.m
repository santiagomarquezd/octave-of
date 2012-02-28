% Solves an 1D Mixture Model using centered schemes for 
% alpha and algebraic momentum equation

page_screen_output(0);

% Initialization

Garcia_Cascales_Phase_separation


dx=xC(2)-xC(1);

% **************************** MAIN PROGRAM ***************************

% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

    % Drift velocity calculation
    % calcVdrp
    calcVdrp

    % U calculation from alphag
    UEqnSimple

    Vg=assign(U,Vdrp,'*');
  
    % Limiting
    % Sweby's fuction calculation
    phiAlphag=superbee(rvalue(alphag,1E-9)); 
    % Limited values calculation
    [alphagLimited]=limitedValues(alphag,phiAlphag,dx,dt);
    % Local velocities at interfaces
    aeigens=simpleSpeed(V0, alphag,rhol, rhog, rhom, N);
    % Time advancement by Kurnanov & Tadmor's scheme
    [alphagDummy,alphag]=KT(alphagh,alphag,alphagLimited,alphagLimited,@alphaSimpleFlux,@alphaSimpleFlux,aeigens,dx,dt,0,rhol,rhog,1);
  end
  
  % Mixture density actualization
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
 
  % BC's actualization
  U=setBC(U,constField(0,N),xC,xF,g);
  alphag=setBC(alphag,rhom,xC,xF,g);
  %rhom=setBC(rhom,constField(0,N),xC,xF,g);
  
end % Temporal loop


