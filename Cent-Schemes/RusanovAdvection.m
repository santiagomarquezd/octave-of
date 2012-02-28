% Solves a simple advection problem by Rosunov's Method
% du/dt+d/dx[F(u)]=0; F(u)=v*u

% Variables clearance
clear all

% Domain extension
xleft=0;
xright=1;

% Physical paramaters
% Advective velocity
vvalue=1;
uvalue=1;

% BC's
uleft=uvalue;
uright=0;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=100;

% Number of cells
N=100;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;

% Eigenvalues
a=vvalue;

% Fields initialization
v=ones(N,1)*vvalue;
u=ones(N,1)*uvalue;
u(ceil(N/2):N)=0;

% Temporal loop
for i=1:timesteps

  % Flux evaluation
  flux=v.*u;
  fluxleft=vvalue*uleft;
  fluxright=vvalue*uright;
  
  % Rosunov scheme
  u(2:end-1)=u(2:end-1)-lambda/2*(flux(3:end)-flux(1:end-2))+1/2*lambda*a*(u(3:end)-u(2:end-1)-u(2:end-1)+u(1:end-2));
   
  % BC's
  u(1)=u(1)-lambda/2*(flux(2)-fluxleft)+1/2*lambda*a*(u(2)-u(1)-u(1)+uleft);
  u(end)=u(end)-lambda/2*(fluxright-flux(end-1))+1/2*lambda*a*(uright-u(end)-u(end)+u(end-1));;

end
