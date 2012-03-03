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
timesteps=2; %100;

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

  if 1  
    % Continuous version (mean velocity calculated by the same formula in all of domain extension)
    
    % Rusanov scheme
    % u_j^(n+1)=u_j^n-lambda/2*(f(u_(j+1)^n)-f(u_(j-1)^n))+1/2*(lambda*a_(j+1/2)^n*(u_(j+1)^n-u_j^n)-lambda*a_(j-1/2)^n*(u_j^n-u_(j-1)^n))   

    % Flux evaluation
    rhomLeft=rhog*alphagLeft+(1-alphagLeft)*rhol;
    fluxleft=(V0*(1-alphagLeft)*(1-alphagLeft*rhog/rhomLeft)+V0*(1-alphagLeft)*alphagLeft.*(rhog/rhomLeft-1))*alphagLeft;
    rhomRight=rhog*alphagRight+(1-alphagRight)*rhol;
    fluxright=(V0*(1-alphagRight)*(1-alphagRight*rhog/rhomRight)+V0*(1-alphagRight)*alphagRight.*(rhog/rhomRight-1)).*alphagRight;
    rhom=rhog*alphag+(1-alphag)*rhol;
    flux=(V0.*(1-alphag(1:end)).*(1-alphag(1:end).*rhog./rhom)+V0.*(1-alphag(1:end)).*alphag(1:end).*(rhog./rhom-1)).*alphag(1:end);
    
    % Local velocities calculation
    a=0.282*(2*alphag-2);
    
    % Local velocities at boundaries
    aLeft=0.282*(2*alphagLeft-2);
    aRight=0.282*(2*alphagRight-2);
    
    a_j_minus_half=[aLeft;max(a(1:end-1),a(2:end))];
    a_j_plus_half=[max(a(1:end-1),a(2:end));aRight];
  
  else
  
    % Discontinuous version (zero mean velocity at alphag=0 and alphag=1 zones)
  
  end
  
  % Rosunov scheme for non boundary cells
  alphag(2:end-1)=alphag(2:end-1)-lambda/2*(flux(3:end)-flux(1:end-2))+1/2*lambda*(a_j_plus_half(2:end-1).*(alphag(3:end)-alphag(2:end-1))-a_j_minus_half(2:end-1).*(alphag(2:end-1)-alphag(1:end-2)));
   
  % BC's
  alphag(1)=alphag(1)-lambda*(flux(2)-fluxleft)+lambda*(a_j_plus_half(1)*(alphag(2)-alphag(1))-a_j_minus_half(1)*(alphag(1)-alphagLeft));
  alphag(end)=alphag(end)-lambda*(fluxright-flux(end-1))+lambda*(a_j_plus_half(end)*(alphagRight-alphag(end))-a_j_minus_half(end)*(alphag(end)-alphag(end-1)));

end
