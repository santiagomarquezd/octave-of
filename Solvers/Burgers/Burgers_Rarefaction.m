% Solves the Burgers equation by Rusanov's Method
% Rarefaction wave case
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[1/2*u^2]
%
% Section 2.4.2 from Riemann Solvers and Numerical Methods for Fluid Dynamics (E.F. Toro)

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
% Dummy value, only for methods compatibility
g=0;

% Domain extension
xleft=-5;
xright=5;

% BC's
uLeft=0;
uRight=1;

% Time-step
dt=0.0001;

% Number of timesteps
timesteps=20000; %100;

% Number of cells
N=10000;

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
  u.internal=ones(N,1)*uRight;
  u.left.type='G';
  u.left.gradient=0;
  u.right.type='G';
  u.right.gradient=0;
  u.internal(1:floor(N/2)+1)=uLeft;
  u=setBC(u,constField(0,N),xC,xF,0);
end

% Fluxes initialization
fluxU.internal=zeros(N,1);
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;

dummyRho=constField(1,N);

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  i

      fluxU=assign(constField(1/2,N),assign(u,constField(2,N),'^'),'*');
      fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

      cfluxU=zeros(N+1,1);
  
      [a_j_minus_half,a_j_plus_half]=aspeedBurgers(u);
  
    
      if 1
	[u,u]=Rusanov(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	   
      else
	[u,u]=LxF(u,u,fluxU,fluxU,dx,dt);
      end

end
