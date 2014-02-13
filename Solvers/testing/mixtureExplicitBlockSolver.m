% Solves an 1D Mixture Model problem using Riemannn-free solvers for the
% U-alpha block

%function mixtureSolver

clear all;
page_screen_output(0);

  Gastaldo_Phase_separation
  %Phase_separation
  %Plume	  

% **************************** MAIN PROGRAM ***************************

% Array of matrices for FVS methods
Aarray=zeros(2,2,N);

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

% Generates Alpha and AlphagRhom auxiliary fields
if 1
  % Creates Alphag0 and Alphag as a copy of alpha
  Alphag0=alphag;
  Alphag=Alphag0;

  % Determination of Alpha auxiliary field
  Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
  Alphag0=setBC(Alphag0,rhom0,xC,xF,g);

  % Creates AlphagRhom0 and AlphagRhom as a copy of alpha.
  AlphagRhom0=Alphag;
  AlphagRhom=AlphagRhom0;

  % Determination of AlphagRhom auxiliary field
  AlphagRhom0.internal=Alphag0.internal.*rhom0.internal;
  AlphagRhom0=setBC(AlphagRhom0,rhom0,xC,xF,g);
end


if 1 % Enables temporal loop

% Temporal loop
for step=1:timesteps
  fprintf('**************** Starting timestep: %d ****************\n',step)

  disp('Masa al comenzar el paso de tiempo (rhom^n)')
  sum(rhom0.internal)

  % Mixture density equation solution
  if (step<TS)
    rhoEqn
  end

  disp('Masa luego de predecir la rhom con la ecuacion de conservacion (rhom~)')
  sum(rhom.internal)
 
  % Drift velocity calculation
  calcVdrp

  %keyboard; pause;

  if 0
    % UEqn and alphagEqn block solving via Rusanov scheme
    % Temporal advancement via Rusanov scheme for U (momentum predictor) and alphag
    if 1
      % Alphag advection like Ishii
      AlphagUBlock
    else
      % alphag and rhom constitutive equation version
      alphagUBlock
    end
  elseif 1
    % Probar haciendo AlphaUBlock con Godunov y usando downwind arriba y centrado abajo a pata
    % Probar el caso de Gastaldo con Vr=cte.

    % Usual calculus for U, Riemmann free for Alpha/Godunov (Flux Vector Splitting) for Alpha

    %disp('Hola')

    % UEqn
    UEqn
    %UEqnVisco
  
    %keyboard; pause;

    % alphaEqn
    %AlphaRF	% Stabilization with Rusanov-like scheme	
    %AlphaTest % Upwind based in front velocities
    %AlphaFVS
    if (step<TS)
      AlphaComp
      %AlphaComp2
    else
      AlphaDelayed
    end
   

    disp('Masa luego de resolver la ecuacion de A (rhom^n+1)')
    sum(rhom.internal)

    %UEqn

  else
    % Traditional mixtureSolver version
  
    % UEqn
    UEqn

    % alphaEqn
    alphaEqnIshii
  
  end

  %PISO loop
  if 1
    for corr=1:nCorr
      if 0
	rhomPhiStill=rhomPhi;
	pEqnExplicitBlock
      elseif 1
	% Traditional PISO
	%alphaEqnIshii
	%AlphaComp2
	pEqn
      elseif 0
	% Vpq actualization
	Vpq=assign(constField(V0,N),assign(constField(1,N),alphag,'-'),'*');
	% Analitycal U
	U=assign(Vpq,assign(alphag,assign(assign(constField(rhog,N),rhom,'/'),constField(1,N),'-'),'*'),'*');
	% BC's setting
	U=setBC(U,constField(0,N),xC,xF,0);
	% Flux updating
	%rhomPhi=fvc_interpolate(assign(rhom,U,'*'), w, xC, xF)
	rhomPhi=fvc_interpolate(rhom, w, xC, xF).*fvc_interpolate(U, w, xC, xF);
      end
    end
  end

  %alphaEqnIshii
  %AlphaComp2

  disp('Masa luego del PISO (rhom^n+1 final)')
  sum(rhom.internal)
  disp('Diferencia de masa')
  deltaRho=sum(rhom.internal)-sum(rhom0.internal)
  disp('Flujo necesario en el borde superior para cerrar la conservacion de masa deltaRho*dx/dt')
  FTeorico=-deltaRho*dx/dt
  disp('Flujo en el borde superior (calculado por PISO)')
  FPISO=rhomPhi(end)
  disp('Error en flujos')
  (FPISO-FTeorico)/(max(FTeorico,FPISO))
  % Temporal variable storing
  TAux(step,1)=sum(rhom.internal); %deltaRho;
%    disp('Velocidad necesaria en el borde superior para cerrar la conservacion de masa rhomPhi/rhom')
%    FTeorico/rhom.right.setvalue
%    disp('Velocidad en el borde superior (calculado por PISO)')
%    U.right.setvalue

if (step<timesteps)

  %keyboard; pause;

  % Set fields as 'old' states
  rhom0bak=rhom0; % For time derivative processing and time shifting in rhom
  U0bak=U0; % For time derivative processing
  rhom0=rhom;
  alphag0=alphag;
  U0=U;
  rhomPhi0=rhomPhi;
  Alphag0=Alphag;
  AlphagRhom0=AlphagRhom;
end

end

end % Ends condition for temporal loop

save -binary dump.dat rhom0 rhom alphag0 alphag U0 U rhomPhi0 rhomPhi Alphag Alphag0 rhom0bak

%end % End function mixtureSolver

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


