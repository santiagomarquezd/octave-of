% Solves two waves of non-linear advection
% using Flux Vector Splitting
% du/dt+d/dx[F(u)]=0; 
% u=[u1;u2], F(u)=[v1*u1;v2*u2]
% Where:
%       vi=vri*(1-ui)-sum_{j!=i}uj*vrj
%       vr=V0i*(1-beta)^ai
%       beta=sum_1^N ui 

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
% Dummy value, only for methods compatibility
g=0;

% Constants for advective velocities
V0=[-0.5;-1;-2];

% Exponents for advective velocities
a=[1;1;1];

% Domain extension
xleft=0;
xright=1;

% BC's
vLeft1=1;
vLeft2=2;

% Initial values
% Two section iniatilization
layers=1;
layerL1=1;
ULeft1=0.4;
layerL2=1;
ULeft2=0.2;
layerL3=1;
ULeft3=0.1;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=2500;

% Number of cells
N=1600;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;

% Fields initialization
% u1
u1.internal=zeros(N,1);
u1.left.type='G';
u1.left.gradient=0;
u1.right.type='G';
u1.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL1);
  u1.internal(1:nCellsLL)=ULeft1;
end
u1=setBC(u1,constField(0,N),xC,xF,0);

% u2
u2.internal=zeros(N,1);
u2.left.type='G';
u2.left.gradient=0;
u2.right.type='G';
u2.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL2);
  u2.internal(1:nCellsLL)=ULeft2;
end
u2=setBC(u2,constField(0,N),xC,xF,0);

% u3
u3.internal=zeros(N,1);
u3.left.type='G';
u3.left.gradient=0;
u3.right.type='G';
u3.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL3);
  u3.internal(1:nCellsLL)=ULeft3;
end
u3=setBC(u3,constField(0,N),xC,xF,0);

% One Jacobian per cell
A=zeros(2,2,N);

% Temporal loop
for i=1:timesteps

    % Prints present timestep
    i

    % Temporal data
    u1tmp=u1.internal;
    u1left=u1.left.setvalue;
    u1right=u1.right.setvalue;

    u2tmp=u2.internal;
    u2left=u2.left.setvalue;
    u2right=u2.right.setvalue;

    u3tmp=u3.internal;
    u3left=u3.left.setvalue;
    u3right=u3.right.setvalue;

    usum=u1.internal+u2.internal+u3.internal;
%      usumleft=u1.left.setvalue+u2.left.setvalue;
%      usumright=u1.right.setvalue+u2.right.setvalue;

    if (0)
    
        % Decoupled
    
        % u1 temporal advancement
    
        % Non boundary cells
        u1.internal(2:N-1)=1/2*(u1tmp(3:N)+u1tmp(1:N-2))-dt/dx/2*(V0(1,1).*(1-u1tmp(3:N)).*(1-u1tmp(3:N)).^a(1,1).*u1tmp(3:N)-V0(1,1).*(1-u1tmp(1:N-2)).*(1-u1tmp(1:N-2)).^a(1,1).*u1tmp(1:N-2));

        % Boundary cells
        % First cell
        u1.internal(1)=1/2*(u1tmp(2)+u1left)-dt/dx/2*(V0(1,1).*(1-u1tmp(2)).*(1-u1tmp(2))^a(1,1).*u1tmp(2)-0); % Impermeable wall
        % Last cell
        u1.internal(N)=1/2*(u1right+u1tmp(N-1))-dt/dx/2*(0-V0(1,1).*(1-u1tmp(N-1)).*(1-u1tmp(N-1))^a(1,1).*u1tmp(N-1)); % Impermeable wall

        % u2 temporal advancement

        % Non boundary cells
        u2.internal(2:N-1)=1/2*(u2tmp(3:N)+u2tmp(1:N-2))-dt/dx/2*(V0(2,1).*(1-u2tmp(3:N)).*(1-u2tmp(3:N)).^a(2,1).*u2tmp(3:N)-V0(2,1).*(1-u2tmp(1:N-2)).*(1-u2tmp(1:N-2)).^a(2,1).*u2tmp(1:N-2));

        % Boundary cells
        % First cell
        u2.internal(1)=1/2*(u2tmp(2)+u2left)-dt/dx/2*(V0(2,1).*(1-u2tmp(2)).*(1-u2tmp(2))^a(2,1).*u2tmp(2)-0); % Impermeable wall
        % Last cell
        u2.internal(N)=1/2*(u2right+u2tmp(N-1))-dt/dx/2*(0-V0(2,1).*(1-u2tmp(N-1)).*(1-u2tmp(N-1))^a(2,1).*u2tmp(N-1)); % Impermeable wall

    
    else

        % Coupled

        % u1 temporal advancement
    
        % Non boundary cells
        u1.internal(2:N-1)=1/2*(u1tmp(3:N)+u1tmp(1:N-2))-dt/dx/2*(V0(1,1).*(1-u1tmp(3:N)).*(1-usum(3:N)).^a(1,1).*u1tmp(3:N)-V0(1,1).*(1-u1tmp(1:N-2)).*(1-usum(1:N-2)).^a(1,1).*u1tmp(1:N-2));

        % Boundary cells
        % First cell
        u1.internal(1)=1/2*(u1tmp(2)+u1left)-dt/dx/2*(V0(1,1).*(1-u1tmp(2)).*(1-usum(2))^a(1,1).*u1tmp(2)-0); % Impermeable wall
        % Last cell
        u1.internal(N)=1/2*(u1right+u1tmp(N-1))-dt/dx/2*(0-V0(1,1).*(1-u1tmp(N-1)).*(1-usum(N-1))^a(1,1).*u1tmp(N-1)); % Impermeable wall

        % u2 temporal advancement

        % Non boundary cells
        u2.internal(2:N-1)=1/2*(u2tmp(3:N)+u2tmp(1:N-2))-dt/dx/2*(V0(2,1).*(1-u2tmp(3:N)).*(1-usum(3:N)).^a(2,1).*u2tmp(3:N)-V0(2,1).*(1-u2tmp(1:N-2)).*(1-usum(1:N-2)).^a(2,1).*u2tmp(1:N-2));

        % Boundary cells
        % First cell
        u2.internal(1)=1/2*(u2tmp(2)+u2left)-dt/dx/2*(V0(2,1).*(1-u2tmp(2)).*(1-usum(2))^a(2,1).*u2tmp(2)-0); % Impermeable wall
        % Last cell
        u2.internal(N)=1/2*(u2right+u2tmp(N-1))-dt/dx/2*(0-V0(2,1).*(1-u2tmp(N-1)).*(1-usum(N-1))^a(2,1).*u2tmp(N-1)); % Impermeable wall

        % u3 temporal advancement

        % Non boundary cells
        u3.internal(2:N-1)=1/2*(u3tmp(3:N)+u3tmp(1:N-2))-dt/dx/2*(V0(2,1).*(1-u3tmp(3:N)).*(1-usum(3:N)).^a(2,1).*u3tmp(3:N)-V0(2,1).*(1-u3tmp(1:N-2)).*(1-usum(1:N-2)).^a(2,1).*u3tmp(1:N-2));

        % Boundary cells
        % First cell
        u3.internal(1)=1/2*(u3tmp(2)+u3left)-dt/dx/2*(V0(2,1).*(1-u3tmp(2)).*(1-usum(2))^a(2,1).*u3tmp(2)-0); % Impermeable wall
        % Last cell
        u3.internal(N)=1/2*(u3right+u3tmp(N-1))-dt/dx/2*(0-V0(2,1).*(1-u3tmp(N-1)).*(1-usum(N-1))^a(2,1).*u3tmp(N-1)); % Impermeable wall

    end

    % Apply BC's
    u1=setBC(u1,constField(0,N),xC,xF,0);

    u2=setBC(u2,constField(0,N),xC,xF,0);

    u3=setBC(u3,constField(0,N),xC,xF,0);

end

close all; plot(xC,u1.internal,'m'); hold on; plot(xC,u2.internal, 'r'); plot(xC,u3.internal, 'c');plot(xC,u1.internal+u2.internal+u3.internal, 'k')

    % Arrays for all cells (one advection matriz per cell)
%      for j=1:N
%          u=[u1.internal(j);u2.internal(j)];
%          A(:,:,j)=pFluxJacobian(u,V0,a);
%      end
%  
%      [V,VT,LAMBDA]=arrayEig(A);
%  
%      % A matrix splittings
%      Aplus=matrixArrayProd(matrixArrayProd(V,lambdaPlus(LAMBDA)),VT);
%      Aminus=matrixArrayProd(matrixArrayProd(V,lambdaMinus(LAMBDA)),VT);  

    % Prints the actual iteration
    
    %[u1,u2]=FVS(u1,u2,Aplus,Aminus,dx,dt,N,1);

    %u2.internal=u2.internal*0;




