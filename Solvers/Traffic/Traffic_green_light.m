% Solves the Traffic equation by Rusanov's Method
% Right red light case
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[V0*(1-u)*u]


% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
% Dummy value, only for methods compatibility
g=0;
V0=1;

% Domain extension
xleft=0;
xright=1;

% BC's and initial conditions
FLeft=0;
FRight=0;
uInit=0.9;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=100; %100;

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
% h
if 0
  h.internal=ones(N,1);
  h.left.type='V';
  h.left.value=hLeft;
  h.right.type='V';
  h.right.value=hRight;
  h.internal(1:floor(N/2)+1)=hLeft;
  h=setBC(h,constField(0,N),xC,xF,0);
else
  u.internal=ones(N,1)*uInit;
  u.left.type='G';
  u.left.gradient=0;
  u.right.type='G';
  u.right.gradient=0;
  u.internal(1:floor(N/2)+1)=0.2;
  u=setBC(u,constField(0,N),xC,xF,0);
end

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
      if 0
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
      
      if 0
	if (rem(i,100)==0)
	  plot(xC,u.internal,'r*-')
	end
      end
end
