% Solves the alpha equation (isolated from the complete multiphase solver)
% with Um!=0 by different methods
%
% u is more dense phase volume fraction (alpha_l)
% 
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u) = [Um - V0*u^a*u*(1 - u)]
% a = 1
% Since Um is divergence free the value selected as initial condition is
% constant over time and space
 
% Variables clearance
clear all;
%close all;
page_screen_output(0);

if 0
    % Perfect sedimentation (convex flux)
    % without layers
    
    % Physical paramaters
    g = 0;
    V0 = 1;
    rhol = 1000;
    rhog = 1;
    % Exponent for relative velocity law
    aexp = 6;

    % Domain extension
    xleft = 0;
    xright = 1;
    
    % Center of volume velocity values
    UmValue = 0;

    % Impermeable wall for fluxes  
    impWall = 1;
    
    % BC's and initial conditions
    FLeft = 0;
    FRight = 0;
    
    % Selects initialization with layers
    layers = 0;

    % Initial u of top u in case of layers
    uTop = 0.5;
    UBottom = 1;
    % Bottom layer height
    layerH = 0.8;
    
elseif 0
    % Perfect sedimentation (convex flux)
    % with layers (PhD. Thesis convex case a)
   
    % Physical paramaters
    g = 0;
    V0 = 1;
    rhol = 1000;
    rhog = 1;
    % Exponent for relative velocity law
    aexp = 0;

    % Domain extension
    xleft = 0;
    xright = 1;
    
    % Center of volume velocity values
    UmValue = 0;

    % No impermeable walls
    impWall = 1;
    
    % BC's and initial conditions
    FLeft = 0;
    FRight = 0;
    
    % Selects initialization with layers
    layers = 1;

    % Initial u of top u in case of layers
    uTop = 0.3;
    UBottom = 0.4;
    % Bottom layer height
    layerH = 0.5;
  
elseif 0
    % Perfect sedimentation (convex flux)
    % with layers (PhD. Thesis convex case b)
   
    % Physical paramaters
    g = 0;
    V0 = 1;
    rhol = 1000;
    rhog = 1;
    % Exponent for relative velocity law
    aexp = 0;

    % Domain extension
    xleft = 0;
    xright = 1;
    
    % Center of volume velocity values
    UmValue = 0;

    % No impermeable walls
    impWall = 1;
    
    % BC's and initial conditions
    FLeft = 0;
    FRight = 0;
    
    % Selects initialization with layers
    layers = 1;

    % Initial u of top u in case of layers
    uTop = 0.5;
    UBottom = 0.3;
    % Bottom layer height
    layerH = 0.5;

elseif 1
    % Non-perfect sedimentation (non-convex flux)
    % (PhD. Thesis non-convex case)
   
    % Physical paramaters
    g = 0;
    V0 = 1;
    rhol = 1000;
    rhog = 1;
    % Exponent for relative velocity law
    aexp = 1;

    % Domain extension
    xleft = 0;
    xright = 1;
    
    % Center of volume velocity values
    UmValue = 0;

    % No impermeable walls
    impWall = 1;
    
    % BC's and initial conditions
    FLeft = 0;
    FRight = 0;
    
    % Selects initialization with layers
    layers = 0;

    % Initial u of top u in case of layers
    uTop = 0.7;
    UBottom = 0.3;
    % Bottom layer height
    layerH = 0.5;    
    
end

% Time-step
dt = 0.001; %0.001/2; %0.001;

% Initialization
initia = 1;

% Method selection, KTcFlux, Godunov
method = 'KTcFlux';

% Inclusion of Um in total flux for Rusanov method
UmIncluded = 1; %1: included, 0: not included

% Number of timesteps
timesteps = 500; %100;

% Number of cells
N = 100; %1000;

% Numerical Pre-processing

% Auxiliar variable for temporal debugging
TAux = zeros(timesteps, 1);

% Cell centers and face centers calculation
% Equi-spaced cells
dx = (xright - xleft)/N;
xC = ((xleft + dx/2):dx:(xright - dx/2))';
xF = (xleft:dx:xright)';

% Interpolation weights calculation
w = weights(xC, xF);

% Face areas
S = 1;
Sf = ones(size(xF))*S;

% Fields initialization
% u
u.internal = ones(N,1)*uTop;
u.left.type = 'G';
u.left.gradient = 0;
u.right.type = 'G';
u.right.gradient = 0;

if 1
  % Particular field initializations
  if 0 
    u.internal(1:1) = 1;
    u.internal(end:end) = 0;
  else
    u.internal = ones(N,1)*uTop;
    u.left.type = 'V';
    u.left.value = 0;
    u.right.type = 'G';
    u.right.gradient=0;
  end
end

if (layers)
    % Detect number of cells in bottom layer
    nCellsBL = sum(xC < layerH);
    u.internal(1:nCellsBL) = UBottom;
end
u = setBC(u, constField(0, N), xC, xF, 0);

% Um
if 1
    Um.internal = zeros(N,1);
    Um.left.type = 'G';
    Um.left.gradient = 0;
    Um.right.type = 'G';
    Um.right.gradient = 0;
    Um = setBC(Um,constField(0,N), xC, xF, 0);
else
    Um.internal = ones(N, 1)*UmValue;
    Um.left.type = 'V';
    Um.left.value = UmValue;
    Um.right.type = 'V';
    Um.right.value = UmValue;
    Um = setBC(Um, constField(0, N), xC, xF, 0);
end

% Fluxes initialization
fluxU.internal = zeros(N, 1);
fluxU.left.type = 'V';
fluxU.left.value = 0;
fluxU.right.type = 'V';
fluxU.right.value = 0;
fluxU = setBC(fluxU, constField(0, N), xC, xF, 0);

% Dummy rho for set BC operations
dummyRho = constField(1, N);
  
% Memory allocation only needed for FVS
if (strcmp(method, 'FVS'))
  A = zeros(2, 2 ,N);
end

% Initialization or not
if (initia != 1)
  disp('Continued running');
  load('data.dat');
end

% Temporal loop
for i = 1:timesteps

    % Prints the actual iteration
    printf('Time-step: %d. Time: %g\n',i,i*dt);

    % Sum of u along the domain to check conservation
    acc = sum(u.internal);	
    printf('Sum of u in the domain: %g\n', acc);

    % Common fields calculation
    Urlg = assign(constField(-V0, N), assign(u, constField(aexp, N), '^'), '*');

    if (strcmp(method, 'Godunov'))

    % No conservative fluxes are used
    cfluxU = zeros(N+1,1);

    [u] = Godunov(u, @alphaEqnUmFlux, cfluxU, dummyRho, dummyRho,...
            V0, rhol, rhog, 1, dx, dt, 10);

    elseif (strcmp(method,'KTcFlux'))

    % Using Kurganov & Tadmor but with Um as centered flux

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag = superbee(rvalue(u,1E-9));
    %phiAlphag = vanLeer(rvalue(u,1E-9));
    % Constant values by cells (mimiking Rusanov?)
    phiAlphag = vanLeer(rvalue(u,1E-9))*0;

    %keyboard; pause;

    % Limited values calculation
    [uLimited] = limitedValues(u, phiAlphag, dx, dt);

    % Eigenvalues for complete flux
    [a_j_minus_half, a_j_plus_half] = ...
        aspeedIsolatedAlphaEqnUm(u, Um, rhol, rhog, V0, aexp);    

    %keyboard; pause;

    % Stabilization flux has to be zero at boundaries
    a_j_minus_half(1) = 0;
    a_j_plus_half(end) = 0;
    % Flatten for KT function
    aeigens = [a_j_minus_half; a_j_plus_half(end)];

    % Um as a centered flux
    phic = fvc_interpolate(Um, w, xC, xF);

    % Impermeable walls
    if impWall
        phic(1) = 0;
        phic(end) = 0;
    end

    % alphag at faces
    uInt = fvc_interpolate(u, w, xC, xF);

    % Time advancement by Kurnanov & Tadmor's scheme
    % with additional centered flux
    [u] = KTcFlux(u, uLimited, @alphaEqnUmFluxFlat, aeigens, phic, uInt,...
            dx, dt, constField(1, N), constField(1, N),...
            ones(N + 1, 1), V0, rhol, rhog, aexp);

    end

    % BC adjusting
    u = setBC(u, constField(0, N), xC, xF, g);

    % Printing and saving
    if 0
        hold on;
        if (rem(i,100) == 0 || i == 1)
            plot(xC,u.internal,'r*-')
            if 0
                eval(['save alphaEqnIsolatedUm-' method '-' num2str(i)...
                '.dat u Um Um0 Urlg Sf xC xF w'])
            end
        end
    end

    % Time debugging variable storing
    TAux(i,1) = sum(u.internal)-acc;
    printf('Delta u in present timestep: %g\n', TAux(i));

    % Field actualization
    Um0 = Um;
    % rhom0=rhom;

end

save data.dat u 

