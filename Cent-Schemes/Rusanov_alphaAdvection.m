% Solves alpha advection equation problem by Rosunov's Method
% du/dt+d/dx[F(u)]=0; F(u)=(V0*(1-alphag)*(1-alphag*rhog/rhom)+
% V0*(1-alphag)*alphag*(rhog/rhom-1))*alphag

% Variables clearance
clear all

% Physical paramaters
rhog=1;
rhol=1000;
V0=0.282;
alphag0=0.5;

% Domain extension
xleft=0;
xright=1;


% BC's
alphagLeft=0;
alphagRight=1;

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

% Fields initialization
alphag=ones(N,1)*alphag0;

% Temporal loop
for i=1:timesteps

  % Flux evaluation
  rhom=rhog*alphag+(1-alphag)*rhol;
  flux=(V0.*(1-alphag(1:end)).*(1-alphag(1:end).*rhog/rhom)+V0.*(1-alphag(1:end)).*alphag(1:end).*(rhog/rhom-1))*alphag(1:end);
  rhomLeft=rhog*alphagLeft+(1-alphagLeft)*rhol;
  fluxleft=(V0*(1-alphagLeft*(1-alphagLeft*rhog/rhom)+V0*(1-alphagLeft)*alphagLeft.*(rhog/rhom-1))*alphagLeft;
  rhomRight=rhog*alphagRight+(1-alphagRight)*rhol;
  fluxright=(V0*(1-alphagRight*(1-alphagRight*rhog/rhomRight)+V0*(1-alphagRight)*alphagRight.*(rhog/rhomRight-1))*alphagRight;
  
  % Local velocities calculation
  
  
  % Rosunov scheme
  alphag(2:end-1)=alphag(2:end-1)-lambda/2*(flux(3:end)-flux(1:end-2))+1/2*lambda*((alphag(3:end)-alphag(2:end-1))-(alphag(2:end-1)-alphag(1:end-2)));
   
  % BC's
  alphag(1)=alphag(1)-lambda/2*(flux(2)-fluxleft)+1/2*lambda*a*(alphag(2)-alphag(1)-alphag(1)+alphagLeft);
  alphag(end)=alphag(end)-lambda/2*(fluxright-flux(end-1))+1/2*lambda*a*(alphagRight-alphag(end)-alphag(end)+alphag(end-1));

end
