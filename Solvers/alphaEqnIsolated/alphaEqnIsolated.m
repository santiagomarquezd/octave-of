% Solves the Alpha equation (isolated from the complete multiphase solver) with Vm=0
% by Rusanov, Flux Vector Splitting UADE and Central Difference Methods
%
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[u*(1-cp)*V0*(1-u)^a], cp=u*rhod/rhom
% rhom=rhod*u+(1-u)*rhol
% a=1

% Variables clearance
clear all;
%close all;
page_screen_output(0);

if 1
  % General case
  % Physical paramaters
  g=0;
  V0=1;
  rhol=1000;
  rhog=1;
  aexp=1;

  % Domain extension
  xleft=0;
  xright=1;
else
  % Latsa/Youngs
  % Physical paramaters
  g=0;
  V0=1;
  rhol=1;
  rhog=0.999;

  % Domain extension
  xleft=0;
  xright=2;
end

% BC's and initial conditions
FLeft=0;
FRight=0;
% Initial u of top u in case of layers
uTop=0.7;%0.5;
layers=1; % Selects initialization with layers
UBottom=0.6;

% Time-step
dt=0.001;

% Initialization
initia=1;

% Method selection, Rusanov, FVS, Centered, Godunov
method='KT';

% Number of timesteps
timesteps=200; %100;

% Number of cells
N=1000/2;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
% Interpolation weights calculation
w=weights(xC, xF);


% Fields initialization
% u
u.internal=ones(N,1)*uTop;
u.left.type='G';
u.left.gradient=0;
u.right.type='G';
u.right.gradient=0;
if (layers)
  u.internal(1:floor(N/2)+1)=UBottom;
end
u=setBC(u,constField(0,N),xC,xF,0);

% Fluxes initialization
fluxU.internal=zeros(N,1);
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;
fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

dummyRho=constField(1,N);

% Memoru allocation only needed for FVS
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
  %i
	 
  if (strcmp(method,'Rusanov'))

    
    rhom=assign(assign(constField(rhog,N),u,'*'),assign(assign(constField(1,N),u,'-'),constField(rhol,N),'*'),'+');
    cp=assign(assign(constField(rhog,N),u,'*'),rhom,'/');
    if 1
      fluxU=assign(assign(assign(constField(V0,N),assign(constField(1,N),u,'-'),'*'),u,'*'),assign(constField(1,N),cp,'-'),'*');
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,constField(0,N));      
    else
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
  
    [u,u]=Rusanov(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	   

  elseif (strcmp(method,'Godunov'))

    % No conservative fluxes are used
    cfluxU=zeros(N+1,1);

    [u]=Godunov(u,@alphaEqnNoVmFlux,cfluxU,dummyRho,dummyRho,V0,rhol,rhog,1,dx,dt,10);

  elseif (strcmp(method,'KT'))

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9));

 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,constField(0,N));  
    % Stabilization flux has to be zero at boundaries
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;
    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Time advancement by Kurnanov & Tadmor's scheme
    [uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnNoVmFluxFlat,@alphaEqnNoVmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);

  else

    Vpq=assign(constField(V0,N),assign(constField(1,N),u,'-'),'*');
    rhom=assign(assign(constField(rhog,N),u,'*'),assign(assign(constField(1,N),u,'-'),constField(rhol,N),'*'),'+');
    cp=assign(assign(constField(rhog,N),u,'*'),rhom,'/');

    if (strcmp(method,'Centered'))
      % FOAM original
      phiVdrp=fvc_interpolate(assign(Vpq,assign(constField(1,N),cp,'-'),'*'), w, xC, xF);
    elseif (strcmp(method,'UADE'))
    % UADE stabilization
      directionFlux=fvc_interpolate(Vpq, w, xC, xF);
      directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
      phiVdrp=fvc_general_interpolate(assign(Vpq,assign(constField(1,N),cp,'-'),'*'), xC, xF,1,directionFlux);
    end

    phiAlpha=phiVdrp;

    % Impermeable condition at walls
    phiAlpha(1)=0;
    phiAlpha(end)=0;

    % Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
    directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
    Alphag0Int=fvc_general_interpolate(u, xC, xF,-1,directionFluxAG0);
    %Alphag0Int=fvc_general_interpolate(u, xC, xF,1,directionFluxAG0);

    % Explicit temporal scheme
    u.internal(1:end)=u.internal(1:end)-dt./dx*(phiAlpha(2:end).*Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1));
  

  end

  % Printing 
  if 0
    hold on;
    if (rem(i,100)==0 || i==1)
      plot(xC,u.internal,'r*-')
    end
  end

end

save data.dat u 


%fluxUm=assign(assign(assign(assign(constField(V0,N),assign(constField(1,N),u,'-'),'*'),u,'*'),assign(assign(assign(constField(rhog,N),rhom,'/'),constField(1,N),'-'),u,'*'),'*'),u,'*');