
% Domain extension
xleft = 0;
xright = 1;

% Cross sectional areas
S = 1;

% Physical paramaters
nu = 1;
rhoData = 1;
g = -10; %9.81;

% Numerical parameters
% Time-step
dt = 1/40;

% Number of timesteps
timesteps = 1;

% Momentum predictor activation (1: enable, 0:disabled)
MomPred = 1;

% PISO corrections
nCorr = 100;

% Number of cells
N = 40;

% Full message selection (1: enable, 0:disabled)
fullVerbose = 0;

% Problem initialization (1: initialization, 0: run from dump data)
initia = 1; 

% ********************* NUMERICAL PRE-PROCESSING ***********************

% Cell centers and face centers calculation
% Equi-spaced cells
dx = (xright - xleft)/N;
xC = ((xleft + dx/2):dx:(xright - dx/2))';
xF = (xleft:dx:xright)';

% Face areas
Sf = ones(size(xF))*S;

% Cell volumes
V = ones(size(xC))*S*dx;

% Interpolation weights calculation
w = weights(xC, xF);

% ********************* INITIAL CONDITIONS ***********************
% BC's: V, value; G: gradient; BP: buoyant pressure

% U field 
U.internal = ones(N,1)*0;
U.left.type ='V';
U.left.value = 0;
U.right.type = 'G';
U.right.gradient = 0;
U = setBC(U, constField(0,N), xC, xF, g);

% Uphi field 
Uphi.internal = ones(N,1)*1;
Uphi.left.type ='V';
Uphi.left.value = 1;
Uphi.right.type = 'V';
Uphi.right.value = 1;
Uphi = setBC(Uphi, constField(0,N), xC, xF, g);

% rho field
rho = constField(rhoData, N);
rho = setBC(rho, constField(0, N), xC, xF, g);

p.internal = ones(N,1)*0;
%p.left.type = 'BP';
p.left.type = 'G';
p.left.gradient = g;
p.right.type = 'V';
p.right.value = 0;
p = setBC(p, rho, xC, xF, g);

% Face viscosities
nuf = fvc_interpolate(constField(nu, N), w, xC, xF).*Sf;
