% Domain extension
xleft=0;
xright=1;
% Number of cells
N=50;

% Cross sectional areas
S=0.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transportProperties
transportProperties.phasea.rho = 1000; %[ 1 -3 0 0 0 ]
transportProperties.phasea.nu = 1e-06; % [ 0 2 -1 0 0 ]
transportProperties.phasea.d =  0.00048; %[ 0 1 0 0 0 0 0 ]

transportProperties.phaseb.rho = 1; %[ 1 -3 0 0 0 ]
transportProperties.phaseb.nu = 1e-06; % [ 0 2 -1 0 0 ]
transportProperties.phaseb.d =  0.0001; %[ 0 1 0 0 0 0 0 ]

transportProperties.Cvm = 0; %[ 0 0 0 0 0 ]
transportProperties.Cl = 0; %[ 0 0 0 0 0 ]
transportProperties.Ct = 0; %[ 0 0 0 0 0 ]
transportProperties.alphaAlpha = 0; %[ 0 0 0 0 0 ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% InterfacialProperties
InterfacialProperties.dragModela = 'SchillerNaumann';
InterfacialProperties.dragModelb = 'SchillerNaumann';
% Drag Models: GidaspowSchillerNaumann
%               SchillerNaumann
%               Ergun
%               GidaspowErgunWenYu
%               SyamlalOBrien
%               WenYu
%               Gibilaro 
InterfacialProperties.dragPhase = 'a';
% dragPhase: a
%            b
%            blended

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dp=0.005*0; % Ensures no relative velocity
g=-10*0;

% Initial values
alpha0 = 0.0;
beta0 = 1;

% Time-step
dt=0.01;

% Number of timesteps
timesteps=1000;

% PISO corrections
nCorr=3;

% Non-orthogonal corrections
nNonOrthCorr=1;




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

% U liquid field 
Ua.internal=ones(N,1);
Ua.left.type='V';
Ua.left.value=1;
Ua.right.type='G';
Ua.right.gradient=0;
Ua=setBC(Ua,constField(0,N),xC,xF,g);

% U gas field 
Ub.internal=ones(N,1);
Ub.left.type='V';
Ub.left.value=1;
Ub.right.type='G';
Ub.right.gradient=0;
Ub=setBC(Ub,constField(0,N),xC,xF,g);


% alpha field 
alpha.internal=ones(N,1)*alpha0;
alpha.left.type='V';
alpha.left.value=1;
alpha.right.type='G';
alpha.right.gradient=0;
alpha=setBC(alpha,constField(0,N),xC,xF,g);

% beta field 
beta = assign(constField(1,N),alpha,'-');

% Set fields as 'old' states
alpha0=alpha;
beta0=beta;
Ua0=Ua;
Ub0=Ub;
DDtUa = zeros(N,1);
DDtUb = zeros(N,1);

% U mixture field 
U = assign(assign(alpha,Ua,'*'),assign(beta,Ub,'*'),'+');

% Phi field initialization from U  fields
phia =fvc_interpolate(Ua0, w, xC, xF).*Sf;
phib =fvc_interpolate(Ub0, w, xC, xF).*Sf;
phia0 = phia;
phib0 = phib;

phi = phia.*fvc_interpolate(alpha, w, xC, xF) + phib.*fvc_interpolate(beta, w, xC, xF);


% Densidad media:
rho = assign(assign(alpha,constField(transportProperties.phasea.rho,N),'*'),assign(beta,constField(transportProperties.phaseb.rho,N),'*'),'+');

% p_rgh field 
p_rgh.internal=ones(N,1)*0;
p_rgh.left.type='G';
p_rgh.left.gradient=0;
p_rgh.right.type='V';
p_rgh.right.value=0;
p_rgh=setBC(p_rgh,rho,xC,xF,g);
