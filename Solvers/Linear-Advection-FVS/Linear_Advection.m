% Solves two simple waves of linear advection
% using Flux Vector Splitting
% du/dt+d/dx[F(u)]=0; 
% u=[u1;u2], F(u)=[a*u1;b*u2]

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
% Dummy value, only for methods compatibility
g=0;

% Advective velocities
a=1;
b=2;

% Domain extension
xleft=0;
xright=10;

% BC's
vLeft1=1;
vLeft2=2;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=1000;

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
% u1
u1.internal=zeros(N,1);
u1.left.type='V';
u1.left.value=vLeft1;
u1.right.type='G';
u1.right.gradient=0;
u1=setBC(u1,constField(0,N),xC,xF,0);

% u2
u2.internal=zeros(N,1);
u2.left.type='V';
u2.left.value=vLeft2;
u2.right.type='G';
u2.right.gradient=0;
u2=setBC(u2,constField(0,N),xC,xF,0);


% Jabobian and splittings are constant in time
Ai=[a 0;0 b];
A=zeros(2,2,N);
% Arrays for all cells
for i=1:N
  A(:,:,i)=Ai;
end
[V,VT,LAMBDA]=arrayEig(A);

% A matrix spllittings
Aplus=matrixArrayProd(matrixArrayProd(V,lambdaPlus(LAMBDA)),VT);
Aminus=matrixArrayProd(matrixArrayProd(V,lambdaMinus(LAMBDA)),VT);

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  i
  
  [u1,u2]=FVS(u1,u2,Aplus,Aminus,dx,dt,N);

end
