% Solves the Traffic equation by Rusanov's Method
% Right red light case
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[V0*(1-u)*u]


% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
V0=1;
% Dummy value, only for methods compatibility
g=0;
rhol=1000;
rhog=1;
a=1;

% Domain extension
xleft=0;
xright=1;

% BC's and initial conditions
FLeft=0;
FRight=0;
uInit=0.7;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=1000/10; %100;

% Number of cells
N=1000;

% Method selection, Rusanov, FVS, Centered
method='UADE';%'SimpleGodunov';%'Rusanov'; %'UADE';

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;
% Interpolation weights calculation
w=weights(xC, xF);

% Fields initialization
u.internal=ones(N,1)*uInit;
u.left.type='G';
u.left.gradient=0;
u.right.type='G';
u.right.gradient=0;
u.internal(1:floor(N/2)+1)=0.6;%uLeft;
u=setBC(u,constField(0,N),xC,xF,0);


% Fluxes initialization
fluxU.internal=zeros(N,1);
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;
fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

dummyRho=constField(1,N);

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  printf('Time-step: %d. Time: %g\n',i,i*dt);

  if (strcmp(method,'Rusanov'))

    if 1
      % Typical traffic flow flux function
      fluxU=assign(assign(constField(V0,N),assign(constField(1,N),u,'-'),'*'),u,'*');
      [a_j_minus_half,a_j_plus_half]=aspeedTraffic(u,V0);
    else
      % Cubic traffic flow flux function
      fluxU=assign(assign(assign(constField(V0,N),assign(constField(1,N),u,'-'),'*'),u,'*'),assign(constField(1,N),u,'-'),'*');
      [a_j_minus_half,a_j_plus_half]=aspeedTrafficCubic(u,V0);
    end

    % Ensure no fluxes at boundaries
    fluxU.left.type='V';
    fluxU.left.value=0;
    fluxU.right.type='V';
    fluxU.right.value=0;
    fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

    % Stabilization flux has to be zero at boundaries too
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;  

    % No conservative fluxes are used
    cfluxU=zeros(N+1,1);

    if 1
      [u,u]=Rusanov(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	   
    else
      [u,u]=LxF(u,u,fluxU,fluxU,dx,dt);
    end

  elseif (strcmp(method,'UADE'))

    % Traffic velocity
    if 0    
      Vpq=assign(constField(V0,N),assign(constField(1,N),u,'-'),'*');
    else
      Vpq=assign(constField(V0,N),assign(assign(constField(1,N),u,'-'),constField(1/2,N),'^'),'*');
    end
    % UADE stabilization
    directionFlux=fvc_interpolate(Vpq, w, xC, xF);
    directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
    phiVdrp=fvc_general_interpolate(Vpq, xC, xF,1,directionFlux);
    
    phiAlpha=phiVdrp;

    % Impermeable condition at walls
    phiAlpha(1)=0;
    phiAlpha(end)=0;

    directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);

    if 1
      % Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
      Alphag0Int=fvc_general_interpolate(u, xC, xF,-1,directionFluxAG0);
    else
      % Alphag0 values at interfaces with full downwind (direction given by phiAlpha)
      Alphag0Int=fvc_general_interpolate(u, xC, xF,1,directionFluxAG0);
    end

    % Explicit temporal scheme
    u.internal(1:end)=u.internal(1:end)-dt./dx*(phiAlpha(2:end).*Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1));

  elseif (strcmp(method,'SimpleGodunov'))

    % Traffic velocity
    Vpq=assign(constField(V0,N),assign(constField(1,N),u,'-'),'*');
	
    directionUp=(fvc_interpolate(u, w, xC, xF)<=0.5)*(2)-1;
    directionUp=directionUp(2:end-1,1);
    F=fvc_general_interpolate(assign(Vpq,u,'*'), xC, xF,-1,directionUp);

    % Impermeable condition at walls
    F(1)=0;
    F(end)=0;
    
    % Explicit temporal scheme
    u.internal(1:end)=u.internal(1:end)-dt./dx*(F(2:end)-F(1:end-1));
  
  elseif (strcmp(method,'Godunov'))

    % No conservative fluxes are used
    cfluxU=zeros(N+1,1);

    [u]=Godunov(u,@quadraticTrafficFlux,cfluxU,dummyRho,dummyRho,V0,rhol,rhog,1,dx,dt,10);

  end

end
