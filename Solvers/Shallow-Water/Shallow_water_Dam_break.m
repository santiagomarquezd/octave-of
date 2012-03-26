% Solves the Shallow Water equations by Rusanov's Method
% Dam-Break case
% du/dt+d/dx[F(u)]=0; 
% u=[u; hu], F(u)=[uh; hu^2+1/2*g*h^2]
% 
% or
%
% u=[u; hu]=[q1; q2], F(q)=[q2;(q2)^2/q1+1/2*g*q1^2]
%
% Example 13.4 from Finite Volume Methods for Hyperbolic Problems (R. Leveque)

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
g=9.81/10;

% Domain extension
xleft=-5;
xright=5;


% BC's
hLeft=3;
hRight=1;

huLeft=0;
huRight=0;

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
  h.internal=ones(N,1);
  h.left.type='G';
  h.left.gradient=0;
  h.right.type='G';
  h.right.gradient=0;
  h.internal(1:floor(N/2)+1)=hLeft;
  h=setBC(h,constField(0,N),xC,xF,0);
end

% hu
if 0
  hu.internal=ones(N,1);
  hu.left.type='V';
  hu.left.value=huLeft;
  hu.right.type='V';
  hu.right.value=huRight;
  hu=setBC(hu,constField(0,N),xC,xF,0);
else
  hu.internal=zeros(N,1);
  hu.left.type='G';
  hu.left.gradient=0;
  hu.right.type='G';
  hu.right.gradient=0;
  hu=setBC(hu,constField(0,N),xC,xF,0);
end

% Fluxes initialization
fluxH.internal=zeros(N,1);
fluxH.left.type='V';
fluxH.left.value=0;
fluxH.right.type='V';
fluxH.right.value=0;

fluxHU.internal=zeros(N,1);
fluxHU.left.type='V';
fluxHU.left.value=0;
fluxHU.right.type='V';
fluxHU.right.value=0;

dummyRho=constField(1,N);

% Vm and Vdrp allocation
%Vm=flux;
%Vdrp=flux;

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  i
  
  
      fluxH=hu;
      fluxH=setBC(fluxH,constField(0,N),xC,xF,0);


      %(q2)^2/q1+1/2*g*q1^2

      fluxHU=assign(assign(assign(hu,constField(2,N),'^'),h,'/'),assign(assign(constField(1/2,N),constField(g,N),'*'),assign(h,constField(2,N),'^'),'*'),'+');
      fluxHU=setBC(fluxHU,constField(0,N),xC,xF,0);

      cfluxH=zeros(N+1,1);
      cfluxHU=zeros(N+1,1);
  
      [a_j_minus_half,a_j_plus_half]=aspeedShallowWater(h,hu,g);
  
    
      if 1
	[h,hu]=Rusanov(h,hu,fluxH,fluxHU,cfluxH,cfluxHU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	   
      else
	[h,hu]=LxF(h,hu,fluxH,fluxHU,dx,dt);
      end

end

u=assign(hu,h,'/');
eig1=u.internal-sqrt(g*h.internal);
eig2=u.internal+sqrt(g*h.internal);
