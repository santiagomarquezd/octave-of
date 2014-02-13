% Solves the alpha equation (isolated from the complete multiphase solver) with Um!=0
% by the KTcFlux method (and eventually another method for comparison)
%
% u is more dense phase volume fraction (alpha_l)
% 
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[Um-V0*u^a*u*(1-u)]
% a=1
% Since Um is divergence free the value selected as initial condition is constant
% over time and space
 
% Variables clearance
clear all;
%close all;
page_screen_output(0);

% ---------------------------------
% Physical and numerical parameters
% ---------------------------------

% 1D Reactor
% Physical paramaters
g=0;
V0=0.4422*0;
rhol=1000;
rhog=1;
% Exponent for relative velocity law
aexp=1*0;

% Domain extension
xleft=0;
xright=1.5/1.5;

% Center of volume velocity values
UmValue=1;

% Impermeable walls?
impWall=0;

% BC's and initial conditions
% Initial u of top u in case of layers
uTop=0.3*0;%0.5;
layers=0; % Selects initialization with layers
uBottom=1;
% Bottom layer height
layerH=1.05;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=4000/4000*500; %/*0+1; %400; %100;

% Number of cells
N=400/40*40;

% Initialization
initia=1;

% Method selection, KT (+/-), KTcFlux (OK), Rusanov (+/-), Godunov (OK), Centered, UADE
method='KTcFlux';%'KTcFlux';

% Numerical Pre-processing
meshPreprocessing

% Auxiliar variable for temporal debugging
TAux=zeros(timesteps,1);

% Fields initialization
% u
u.internal=ones(N,1)*uTop;
u.left.type='V';
u.left.value=1;
u.right.type='G';
u.right.gradient=0;
if (layers)
  % Detect number of cells in bottom layer
  nCellsBL=sum(xC<layerH);
  u.internal(1:nCellsBL)=uBottom;
end
u=setBC(u,constField(0,N),xC,xF,0);

% Um
% Passing velocity
Um.internal=ones(N,1)*UmValue;
Um.left.type='V';
Um.left.value=UmValue;
Um.right.type='V';
Um.right.value=UmValue;
Um=setBC(Um,constField(0,N),xC,xF,0);

% ----------------------------------------
% End of physical and numerical parameters
% ----------------------------------------

% Dummy rho for set BC operations
dummyRho=constField(1,N);

% Initialization or not
if (initia!=1)
  disp('Continued running');
  load('data.dat');
end

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  printf('Time-step: %d. Time: %g\n',i,i*dt);
  
  % Sum of u along the domain to check conservation
  acc=sum(u.internal);	
  printf('Sum of u in the domain: %g\n',acc);

  % Common fields calculation
  % Relative velocity for bubble model Urlg=-V0*u^a
  Urlg=assign(constField(-V0,N),assign(u,constField(aexp,N),'^'),'*');
	 
  if (strcmp(method,'Rusanov'))

  elseif (strcmp(method,'Godunov'))

  elseif (strcmp(method,'KT'))

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9));
    %phiAlphag=vanLeer(rvalue(u,1E-9))*0; % Constant values by cells (mimiking Rusanov?)
 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    %[a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,constField(0,N));
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnUm(u,Um,rhol,rhog,V0,aexp);    

    % Impermeable walls
    if impWall
      % Stabilization flux zero at boundaries
      a_j_minus_half(1)=0;
      a_j_plus_half(end)=0;
    end

    % Stabilization flux has to be zero at boundaries
    % a_j_minus_half(1)=0;
    % a_j_plus_half(end)=0;

    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Time advancement by Kurnanov & Tadmor's scheme
    %[uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnUmFluxFlat,@alphaEqnUmFluxFlat,Um,aeigens,dx,dt,V0,rhol,rhog,aexp);
    [u,uDummy]=KT(u,u,uLimited,uLimited,@alphaEqnUmFluxFlat,@alphaEqnUmFluxFlat,Um,aeigens,dx,dt,V0,rhol,rhog,aexp);

  elseif (strcmp(method,'KTcFlux'))

    % Using Kurganov & Tadmor but with Um as centered flux

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    % phiAlphag=vanLeer(rvalue(u,1E-9));
    % phiAlphag=minmod(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9))*0; % Constant values by cells (mimiking Rusanov?)

    %keyboard; pause;
 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    % Um as a centered flux
    phic=fvc_interpolate(Um, w, xC, xF);

    if 0 
      % Eigenvalues for complete flux
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnUm(u,Um,rhol,rhog,V0,aexp);    
    elseif 0
      % Eigenvalues for flux without Um
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,aexp);
    elseif 1
      % Eigenvalues for complete flux using phic
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnPhiUm(uLimited,phic,rhol,rhog,V0,aexp);    
    end

    %keyboard; pause;

    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Impermeable walls
    if impWall
      phic(1)=0;
      phic(end)=0;
      % Stabilization flux zero at boundaries
      a_j_minus_half(1)=0;
      a_j_plus_half(end)=0;
    end

    % alphag at faces
    if 0
      % Using linear interpolation is only correct if a constant reconstruction is used. Thesis Eqn. (4.20)
      uInt=fvc_interpolate(u, w, xC, xF);
    else
      % The general implementation implies the mean of the reconstructed values at both sides. Thesis Eqn. (4.24)
      uInt=[uLimited.u_i_m_h_l(1);
	    (uLimited.u_i_p_h_l(1:(end-1))+uLimited.u_i_p_h_r(1:(end-1)))/2;
	    uLimited.u_i_p_h_r(end)];
    end


    % Time advancement by Kurnanov & Tadmor's scheme with additional centered flux
    %[uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnNoUmFluxFlat,@alphaEqnNoUmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);
    [u]=KTcFlux(u,uLimited,@alphaEqnUmFluxFlat,aeigens,phic,uInt,dx,dt,constField(1,N),constField(1,N),ones(N+1,1),V0,rhol,rhog,aexp);

  else
    
   % Dummy option (previous UADE)	  

  end

  % BC adjusting
  u=setBC(u,constField(0,N),xC,xF,g);

  % Printing and saving
  if 0
    hold on;
    if (rem(i,100)==0 || i==1)
      %plot(xC,u.internal,'r*-')
      eval(['save alphaEqnIsolatedUm-' method '-' num2str(i) '.dat u Um Um0 Urlg Sf xC xF w'])
    end
  end

  % Time debugging variable storing
  TAux(i,1)=sum(u.internal)-acc;
  printf('Delta u in present timestep: %g\n',TAux(i));
  
  % Field actualization
  Um0=Um;
  % rhom0=rhom;

end

% Solution backup
save data.dat u 

