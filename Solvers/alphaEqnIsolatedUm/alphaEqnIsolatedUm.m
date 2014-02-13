% Solves the alpha equation (isolated from the complete multiphase solver) with Um!=0
% by different methods
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

if 0
  % General case
  % Physical paramaters
  g=0;
  V0=1;
  rhol=1000;
  rhog=1;
  % Exponent for relative velocity law
  aexp=1;

  % Domain extension
  xleft=0;
  xright=1;
  
  % Center of volume velocity values
  UmValue=0;

  % Impermeable wall for fluxes  
  impWall=1;
elseif 1
  % 1D Reactor
  % Physical paramaters
  g=0;
  V0=1;
  rhol=1000;
  rhog=1;
  % Exponent for relative velocity law
  aexp=1;

  % Domain extension
  xleft=0;
  xright=1;
  
  % Center of volume velocity values
  UmValue=0.4422;

  % No impermeable walls
  impWall=0;
elseif 0
  % Latsa/Youngs
  % Physical parameters
  g=0;
  V0=1;
  rhol=1;
  rhog=0.999;

  % Domain extension
  xleft=0;
  xright=2;

  % Impermeable wall for fluxes  
  impWall=1;
end

% BC's and initial conditions
FLeft=0;
FRight=0;
% Initial u of top u in case of layers
uTop=0;%0.5;
layers=1; % Selects initialization with layers
UBottom=1;
% Bottom layer height
layerH=0.8;

% Time-step
dt=0.001; %0.001/2; %0.001;

% Initialization
initia=1;

% Method selection, KT (+/-), KTcFlux (OK), Rusanov (+/-), Godunov (OK), Centered, UADE
method='KTcFlux';

% Inclusion of Um in total flux for Rusanov method
UmIncluded=1; %1: included, 0: not included

% Number of timesteps
timesteps=50; %400; %100;

% Number of cells
N=100; %1000;

% Numerical Pre-processing

% Auxiliar variable for temporal debugging
TAux=zeros(timesteps,1);

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
% Interpolation weights calculation
w=weights(xC, xF);
% Face areas
S=1;
Sf=ones(size(xF))*S;

% Fields initialization
% u
u.internal=ones(N,1)*uTop;
u.left.type='G';
u.left.gradient=0;
u.right.type='G';
u.right.gradient=0;
if 1
  % Particular field initializations
  if 0 
    u.internal(1:1)=1;
    u.internal(end:end)=0;
  else
    u.internal=ones(N,1)*uTop;
    u.left.type='V';
    u.left.value=0;
    u.right.type='G';
    u.right.gradient=0;
  end
end
if (layers)
  % Detect number of cells in bottom layer
  nCellsBL=sum(xC<layerH);
  u.internal(1:nCellsBL)=UBottom;
end
u=setBC(u,constField(0,N),xC,xF,0);

% Um
if 0
  Um.internal=zeros(N,1);
  Um.left.type='G';
  Um.left.gradient=0;
  Um.right.type='G';
  Um.right.gradient=0;
  Um=setBC(Um,constField(0,N),xC,xF,0);
else
  Um.internal=ones(N,1)*UmValue;
  Um.left.type='V';
  Um.left.value=UmValue;
  Um.right.type='V';
  Um.right.value=UmValue;
  Um=setBC(Um,constField(0,N),xC,xF,0);
end

% Fluxes initialization
fluxU.internal=zeros(N,1);
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;
fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

% Dummy rho for set BC operations
dummyRho=constField(1,N);
  
% Memory allocation only needed for FVS
if (strcmp(method,'FVS'))
  A=zeros(2,2,N);
end

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
  Urlg=assign(constField(-V0,N),assign(u,constField(aexp,N),'^'),'*');

  % Um is calculated from an analytical expression since the complete solver is not available
  % Um is zero in this problem
  % DEPRECATED. Um is constant over time and space in 1D since div(Um)=0.
  % Um=constField(0,N);
  % Um=setBC(Um,constField(0,N),xC,xF,0);
	 
  if (strcmp(method,'Rusanov'))

    if (UmIncluded)
      % Um effect included in total flux
      fluxU=assign(assign(Urlg,assign(u,assign(constField(1,N),u,'-'),'*'),'*'),assign(Um,u,'*'),'+');
      % No conservative/centered fluxes are used
      cfluxU=zeros(N+1,1);
    else
      % Um effect included included as a centered flux in Rusanov scheme
      % In the real solver (Um)_f is given by the PISO loop
      if 1
	% Like FOAM
	cfluxU=fvc_interpolate(Um, w, xC, xF).*fvc_interpolate(u, w, xC, xF);
      else
	cfluxU=fvc_interpolate(assign(Um,u,'*'),w, xC, xF);
      end
      % Zero centered fluxes at boundaries
      fluxU=assign(Urlg,assign(u,assign(constField(1,N),u,'-'),'*'),'*');
    end
    
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnUm(u,rhol,rhog,V0,aexp);  
    % Stabilization flux has to be zero at boundaries too
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;  
  
    % Ensure no fluxes at boundaries
    fluxU.left.type='V';
    fluxU.left.value=0;
    fluxU.right.type='V';
    fluxU.right.value=0;
    fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

      
    % Time advancement
    if (UmIncluded)
      [u,u]=Rusanov(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	
    else
      [u]=MarquezNigro(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,w,xC,xF,dx,dt);
    end

  elseif (strcmp(method,'Godunov'))

    % No conservative fluxes are used
    cfluxU=zeros(N+1,1);

    [u]=Godunov(u,@alphaEqnUmFlux,cfluxU,dummyRho,dummyRho,V0,rhol,rhog,1,dx,dt,10);

  elseif (strcmp(method,'KT'))

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9));
    % phiAlphag=vanLeer(rvalue(u,1E-9))*0; % Constant values by cells (mimiking Rusanov?)
 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    %[a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,constField(0,N));
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnUm(u,Um,rhol,rhog,V0,aexp);    
    % Stabilization flux has to be zero at boundaries
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;
    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Time advancement by Kurnanov & Tadmor's scheme
    [uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnUmFluxFlat,@alphaEqnUmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);

  elseif (strcmp(method,'KTcFlux'))

    % Using Kurganov & Tadmor but with Um as centered flux

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9));
    % phiAlphag=vanLeer(rvalue(u,1E-9))*0; % Constant values by cells (mimiking Rusanov?)

    %keyboard; pause;
 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    if 1  
      % Eigenvalues for complete flux
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnUm(u,Um,rhol,rhog,V0,aexp);    
    else
      % Eigenvalues for flux without Um
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,aexp); 
    end

    %keyboard; pause;

    % Stabilization flux has to be zero at boundaries
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;
    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Um as a centered flux
    phic=fvc_interpolate(Um, w, xC, xF);

    % Impermeable walls
    if impWall
      phic(1)=0;
      phic(end)=0;
    end

    % alphag at faces
    uInt=fvc_interpolate(u, w, xC, xF);

    % Time advancement by Kurnanov & Tadmor's scheme with additional centered flux
    %[uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnNoUmFluxFlat,@alphaEqnNoUmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);
    [u]=KTcFlux(u,uLimited,@alphaEqnUmFluxFlat,aeigens,phic,uInt,dx,dt,constField(1,N),constField(1,N),ones(N+1,1),V0,rhol,rhog,aexp);

  else
    % Centered flux for Um inclusion
    phiAlpha=fvc_interpolate(Um, w, xC, xF);

    if (strcmp(method,'Centered'))
      % FOAM original
      phiVdrp=fvc_interpolate(assign(Urlg,assign(constField(1,N),u,'-'),'*'), w, xC, xF);
    elseif (strcmp(method,'UADE'))
    % UADE stabilization
      directionFlux=fvc_interpolate(Urlg, w, xC, xF);
      directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
      phiVdrp=fvc_general_interpolate(assign(Urlg,assign(constField(1,N),u,'-'),'*'), xC, xF,1,directionFlux);
    end

    phiAlpha=phiVdrp+phiAlpha;

    % Impermeable condition at walls
    phiAlpha(1)=0;
    phiAlpha(end)=0;

    % alphal0 values at interfaces with full upwind (direction given by phiAlpha)
    directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
    alphal0Int=fvc_general_interpolate(u, xC, xF,-1,directionFluxAG0);

    % Explicit temporal scheme
    u.internal(1:end)=u.internal(1:end)-dt./dx*(phiAlpha(2:end).*alphal0Int(2:end)-phiAlpha(1:end-1).*alphal0Int(1:end-1));
  

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

save data.dat u 

