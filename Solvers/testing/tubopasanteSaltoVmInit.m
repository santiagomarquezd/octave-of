clear all

% Domain extension
xleft=0;
xright=1;

% Cross sectional areas
S=0.02;

% Physical paramaters
rhol=10;
mul=0.001;
rhog=1;
dp=0.005*0; % Ensures no relative velocity
g=-10*0;

% Initial values
alphag0=1;

% Time-step
dt=0.0001;

% Number of timesteps
timesteps=0; %40;

% PISO corrections
nCorr=3;

% Non-orthogonal corrections
nNonOrthCorr=0;

% Number of cells
N=100;

% alpha/Vm solution scheme
alphaVmScheme='KT'; %'KT';'LxF'; 'Rusanov'

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
U.internal=10*ones(N,1);
U.internal(51:100)=1;
U.left.type='V';
U.left.value=10;
U.right.type='G';
U.right.gradient=0;
U=setBC(U,constField(0,N),xC,xF,g);

% Vpq field from U field
Vpq=U;
Vpq.left.type='G';
Vpq.left.gradient=0;
Vpq.right.type='G';
Vpq.right.gradient=0;
Vpq=setBC(Vpq,constField(0,N),xC,xF,g);

% alpha field -------------->>>>>>>>>> VER CONDICIONES DE BORDE CORRECTAS
alphag.internal=ones(N,1)*alphag0;
alphag.left.type='V';
alphag.left.value=1;
alphag.right.type='G';
alphag.right.gradient=0;
alphag=setBC(alphag,constField(0,N),xC,xF,g);

% rhom field
rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
rhom=setBC(rhom,constField(0,N),xC,xF,g);

% p_rgh field 
p_rgh.internal=ones(N,1)*0;
p_rgh.left.type='G';
p_rgh.left.gradient=0;
p_rgh.right.type='V';
p_rgh.right.value=0;
p_rgh=setBC(p_rgh,rhom,xC,xF,g);
