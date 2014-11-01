
% Domain extension
xleft=0;
xright=7.5;

% Cross sectional areas
S=1; %0.02;

% Physical paramaters
nul=1E-6;
rhol=1000;
mul=0.001;
nug=1E-5;
rhog=1.2;
dp=0.005; % 
g=-9.81;

% Relative velocity model (1: UADE, 2: Schiller-Naumann, 3:Constant, V0)
VpqModel=3;
% UADE model Vpq=V0.*((alphaMax-min(alpha.internal,alphaMax))/alphaMax).^a
V0=1; %1;%0.282*0;
alphaMax=1;
aexp=0;

% Initial alphag. Top alphag in case of layers
alphaG0Top=0.5;%1;%0.5;
layers=0; % Selects initialization with layers
alphaG0Bottom=0;


% Numerical parameters
% Time-step
dt=0.001/5;%0.141843971631206/10; %0.141843971631206/10; %0.00001;

% Number of timesteps
timesteps=2500; %5000/50; %12;%2000;

% Number of iteration when time-shifting in rhom starts
TS=2; %3E10;

% Numerical diffusivity for stabilization multiplier in Alpha Equation
mult=0; %1    %nu=mult*1/2*mean(abs(U)+abs(Vpq))*mean(dx)*ones(size(rhomPhi));

% Multiplier for Numerical diffusivity in U Equation stabilization
multU=1;

% P-V coupling method 0: PISO, 1:Chorin
PV=0;
% PISO corrections
nCorr=3;

% Non-orthogonal corrections
nNonOrthCorr=1;

% Discretization scheme for alpha divergence term
alphaDiv=1; %2;
% Explicit solution of alpha equation
alphaExplicit=1; %1: explicit, 0: implicit

% Number of cells
N=400; %400; %6;

% Number of PISO correction in the last timestep
stopCorr=3; 

% Full message selection
fullVerbose=0;  %1: enable, 0:disabled

% Problem initialization
initia=0; % 1: initialization, 0: run from dump data

% Temporal auxiliar variable, only for debugging
TAux=zeros(timesteps,1);

% ********************* NUMERICAL PRE-PROCESSING ***********************

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';

% Face areas
Sf=ones(size(xF))*S;

% Cell volumes
V=ones(size(xC))*S*dx;

% Interpolation weights calculation
w=weights(xC, xF);


% ********************* INITIAL CONDITIONS ***********************
% BC's: V, value; G: gradient; BP: bouyant pressure

% U field 
U.internal=0*ones(N,1)*-V0/4;
U.left.type='V';
U.left.value=0;
if 0
  U.right.type='V';
  U.right.value=0;
else
  U.right.type='G';
  U.right.gradient=0;
end
% Detailed initialization
%U.internal(end-3:end)=0;
%U.right.gradient=0;
U=setBC(U,constField(0,N),xC,xF,g);

% Vpq field from U field
Vpq=U;
Vpq.left.type='V';
Vpq.left.value=0;
Vpq.right.type='V';
Vpq.right.value=0;
Vpq=setBC(Vpq,constField(0,N),xC,xF,g);

% alpha field -------------->>>>>>>>>> VER CONDICIONES DE BORDE CORRECTAS
alphag.internal=ones(N,1)*alphaG0Top;
alphag.left.type='V';
alphag.left.value=0;
%alphag.left.gradient=0;
alphag.right.type='G';
%alphag.right.value=1;
alphag.right.gradient=0;
%alphag.internal(1)=0;
%alphag.internal(end)=1;
if (layers)
  alphag.internal(1:floor(N/2)+1)=alphaG0Bottom;
end
%alphag.internal(1:4)=0;
%alphag.internal(end-4:end)=1;
alphag=setBC(alphag,constField(0,N),xC,xF,g);

% rhom field
rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
rhom=setBC(rhom,constField(0,N),xC,xF,g);

if (PV==0)
  % p_rgh field 
  p_rgh.internal=ones(N,1)*0; %5005;
  p_rgh.left.type='BP';
  p_rgh.right.type='V';
  p_rgh.right.value=0; %5005;
  p_rgh=setBC(p_rgh,rhom,xC,xF,g);
elseif (PV==1)
  p.internal=ones(N,1)*0; %5005;
  p.left.type='BP';
  p.right.type='V';
  p.right.value=0;
  p=setBC(p,rhom,xC,xF,g);
end

