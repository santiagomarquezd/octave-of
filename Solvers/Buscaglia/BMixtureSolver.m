% Solves an 1D Mixture Model problem
% Developing with Gustavo C. Buscaglia

clear all;
page_screen_output(0);

  %Garcia_Cascales_Phase_separation
  Phase_separation

% **************************** MAIN PROGRAM ***************************

% Gravity treatment terms
ghf = g*xF;
gh = g*xC;

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

% Temporal loop
for step=1:timesteps
  fprintf('---------------------------------\n')
  fprintf('Timestep: %d. Time: %g\n',step,step*dt)

  % Mixture equation solution
  rhoEqn

  % Drift velocity calculation
  calcVdrp  

if (PV==0)

  % UEqn. U predictor
  UEqn
  
  % alphaEqn
  BAlphaRF

  %if (step<timesteps+1)

    %PISO loop
    if 1  
      for corr=1:nCorr
	%BpEqn
	pEqn
      end
    end

  %end

elseif (PV==1)

  % Chorin's methods
  BUEqnChorin

  %alphaEqnIshii
  BAlphaRFChorin

  BpEqnChorin 

end

  if 1
    % Set fields as 'old' states
    rhom0bak=rhom0; % For time derivative processing
    U0bak=U0; % For time derivative processing
    rhom0=rhom;
    alphag0=alphag;
    U0=U;
    rhomPhi0=rhomPhi;
  end

end

% Save data for restarting
save -binary dump.dat rhom0 rhom alphag0 alphag U0 U rhomPhi0 rhomPhi %Alphag0 Alphag


%  if 0
%    temporal=fvc_ddt(assign(U,rhom,'*'), assign(U0bak,rhom0bak,'*'), dt);
%    convec=fvc_div_face(rhomPhi,V);
%    drift=fvc_div_cell(arg, w, xC, xF, Sf, V);
%    volumetric=fvc_reconstruct((-ghf.*fvc_snGrad(rhom,xC,xF)-fvc_snGrad(p_rgh,xC,xF)).*Sf,Sf);
%  else
%    temporal=ddtM*UmomPred.internal-ddtRHS;
%    convec=convM*UmomPred.internal-convRHS;
%    drift=driftRHS;
%    volumetric=volRHS;
%    residual=temporal+convec-drift-volumetric;
%  end


