% Solves alpha equation for 1D Mixture Model using centered schemes from
% Kurganov and Tadmor with superbee reconstruction

page_screen_output(0);

% Initialization

%Garcia_Cascales_Phase_separation
Phase_separation

% **************************** MAIN PROGRAM ***************************

% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

    % Limiting
    % Sweby's fuction calculation
    phiAlphag=superbee(rvalue(alphag,1E-9)); 
    % Limited values calculation
    [alphagLimited]=limitedValues(alphag,phiAlphag,dx,dt);
    if 0
      % Old version
      % Local velocities at interfaces
      aeigens=simpleSpeed(V0, alphag,rhol, rhog, rhom, N);
    else 
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(alphag,rhol,rhog,V0,constField(0,N));  
      % Stabilization flux has to be zero at boundaries too
      a_j_minus_half(1)=0;
      a_j_plus_half(end)=0;
      % Flatten for KT function
      aeigens=[a_j_minus_half; a_j_plus_half(end)];
    end

    % Time advancement by Kurnanov & Tadmor's scheme
    [alphagDummy,alphag]=KT(alphag,alphag,alphagLimited,alphagLimited,@alphaEqnVmFluxFlat,@alphaEqnVmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);
  %end
  
  alphag=setBC(alphag,rhom,xC,xF,g);
  
end % Temporal loop


