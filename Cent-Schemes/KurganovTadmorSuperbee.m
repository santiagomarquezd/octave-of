% Solves a simple advection problem by Kurganov & Tadmor Scheme
% with Superbee limiter
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

% Variable allocation
r=zeros(N,1);
u_i_plus_half_left=zeros(N,1);
u_i_plus_half_right=zeros(N,1);
u_i_minus_half_left=zeros(N,1);
u_i_minus_half_right=zeros(N,1);

% Temporal loop
for i=1:timesteps

  % r variable calculation for TVD limiter
  % Internal
  r(2:end-1)=(u(2:end-1)-u(1:end-2))./(u(3:end)-u(2:end-1)+1E-9);
  % BC's
  r(1)=2*(u(1)-uleft)/(u(2)-u(1)+1E-9);
  r(end)=(u(end)-u(end-1))/(uright-u(end)+1E-9)/2;
  
  % Limiting
  phi=superbee(r);
  
  % Limited values at interfaces
  u_i_plus_half_left(1:end-1)=u(1:end-1)+0.5*phi(1:end-1).*(u(2:end)-u(1:end-1));
  % Zero gradient BC
  u_i_plus_half_left(end)=u(end);
   
  u_i_plus_half_right(1:end-2)=u(2:end-1)-0.5*phi(2:end-1).*(u(3:end)-u(2:end-1));
  u_i_plus_half_right(end-1)=u(end)-0.5*phi(end)*(uright-u(end));
  % Zero gradient BC
  u_i_plus_half_right(end)=uright;
  
  u_i_minus_half_left(2:end)=u(1:end-1)+0.5*phi(1:end-1).*(u(2:end)-u(1:end-1));
  % Fixed value BC
  u_i_minus_half_left(1)=uleft;
    
  u_i_minus_half_right(1:end-1)=u(1:end-1)-0.5*phi(1:end-1).*(u(2:end)-u(1:end-1));
  % Zero gradient BC
  u_i_minus_half_right(end)=u(end);
  
   
  % Kurganov and Tadmor scheme (simplified for constant scalar advection
  u=u-lambda*vvalue*(u_i_plus_half_left-u_i_minus_half_left);

end
