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
V0=[-2;0];

% Exponents for advective velocities
a=[1;1];

% Domain extension
xleft=0;
xright=1;

% BC's
vLeft1=1;
vLeft2=2;

% Two section iniatilization
layers=1;
layerL1=1;
ULeft1=0.25;
layerL2=1;
ULeft2=0.5;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=1; %1000/1000*10;

% Number of cells
N=4; %400

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
u1.left.type='V';
u1.left.value=vLeft1;
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
u2.left.type='V';
u2.left.value=vLeft2;
u2.right.type='G';
u2.right.gradient=0;
if (layers)
  % Detect number of cells in left zone
  nCellsLL=sum(xC<layerL2);
  u2.internal(1:nCellsLL)=ULeft2;
end
u2=setBC(u2,constField(0,N),xC,xF,0);


% One Jacobian per cell
A=zeros(2,2,N);

% Temporal loop
for i=1:timesteps

    % Arrays for all cells (one advection matriz per cell)
    for j=1:N
        u=[u1.internal(j);u2.internal(j)];
        A(:,:,j)=pFluxJacobian(u,V0,a);
    end

    [V,VT,LAMBDA]=arrayEig(A);

    % A matrix splittings
    Aplus=matrixArrayProd(matrixArrayProd(V,lambdaPlus(LAMBDA)),VT);
    Aminus=matrixArrayProd(matrixArrayProd(V,lambdaMinus(LAMBDA)),VT);  

    % Prints the actual iteration
    i
    
    %[u1,u2]=FVS(u1,u2,Aplus,Aminus,dx,dt,N,1);

    %u2.internal=u2.internal*0;

end
