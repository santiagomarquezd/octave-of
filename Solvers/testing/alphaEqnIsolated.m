% Solves the Alpha equation (isolated from the complete multiphase solver) with Vm=0
% by Rusanov, Flux Vector Splitting and Central Difference Methods
%
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[u*(1-cp)*V0*(1-u)^a], cp=u*rhod/rhom
% rhom=rhod*u+(1-u)*rhol
% a=1

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
g=0;
V0=1;
rhol=1000;
rhog=1;

% Domain extension
xleft=0;
xright=1;

% BC's and initial conditions
FLeft=0;
FRight=0;
uInit=0.5;

% Time-step
dt=0.001;

% Method selection, Rusanov, FVS, Centered
method='Rusanov';

% Number of timesteps
timesteps=1000; %100;

% Number of cells
N=1000;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;

% Fields initialization
% u
u.internal=ones(N,1)*uInit;
u.left.type='G';
u.left.gradient=0;
u.right.type='G';
u.right.gradient=0;
%u.internal(1:floor(N/2)+1)=0.2;
u=setBC(u,constField(0,N),xC,xF,0);

% Fluxes initialization
fluxU.internal=zeros(N,1);
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;
fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

dummyRho=constField(1,N);

hold on;

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  i
  
  if (strcmp(method,'Rusanov'))

    % Cubic traffic flow flux function
    rhom=assign(assign(constField(rhog,N),u,'*'),assign(assign(constField(1,N),u,'-'),constField(rhol,N),'*'),'+');
    cp=assign(constField(1,N),assign(assign(constField(rhog,N),u,'*'),rhom,'/'),'-');
    fluxU=assign(assign(assign(constField(V0,N),assign(constField(1,N),u,'-'),'*'),u,'*'),assign(constField(1,N),cp,'-'),'*');
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,constField(0,N));

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

   elseif (strcmp(method,'FVS'))


   elseif (strcmp(method,'Centered'))

   end


    if 0
      if (rem(i,100)==0)
	plot(xC,u.internal,'r*-')
      end
    end
end
