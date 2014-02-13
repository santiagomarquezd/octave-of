% Solves the Alpha equation (isolated from the complete multiphase solver) with Vm!=0
% by Rusanov, Flux Vector Splitting UADE and Central Difference Methods
%
% du/dt+d/dx[F(u)]=0; 
% u=[u], F(u)=[u*(1-cp)*V0*(1-u)^a], cp=u*rhod/rhom
% rhom=rhod*u+(1-u)*rhol
% a=1

% Variables clearance
clear all;
%close all;
page_screen_output(0);

if 1
  % General case
  % Physical paramaters
  g=0;
  V0=1;
  rhol=1000;
  rhog=1;
  % Exponent for relative velocity law
  aexp=1;

  % Domain extension
  xleft=0;
  xright=1;
else
  % Latsa/Youngs
  % Physical parameters
  g=0;
  V0=1;
  rhol=1;
  rhog=0.999;

  % Domain extension
  xleft=0;
  xright=2;
end

% BC's and initial conditions
FLeft=0;
FRight=0;
% Initial u of top u in case of layers
uTop=0.3;%0.5;
layers=0; % Selects initialization with layers
UBottom=0.9;


% Time-step
dt=0.001/2; %0.001;

% Initialization
initia=1;

% Method selection, KT, KTcFlux, Rusanov, Godunov, Centered, UADE
method='Rusanov';

% Inclusion of Vm in total flux for Rusanov method
VmIncluded=1; %1: included, 0: not included

% Number of timesteps
timesteps=1000; %400; %100;

% Number of cells
N=400; %1000;

% Numerical Pre-processing

% Auxiliar variable for temporal debugging
TAux=zeros(timesteps,1);

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
% Interpolation weights calculation
w=weights(xC, xF);
% Face areas
S=1;
Sf=ones(size(xF))*S;

% Fields initialization
% u
u.internal=ones(N,1)*uTop;
u.left.type='G';
u.left.gradient=0;
u.right.type='G';
u.right.gradient=0;
if (layers)
  u.internal(1:floor(N/2)+1)=UBottom;
end
u.internal(1:4)=0;
u.internal(end-4:end)=1;
u=setBC(u,constField(0,N),xC,xF,0);

% Vm
Vm.internal=zeros(N,1);
Vm.left.type='G';
Vm.left.gradient=0;
Vm.right.type='G';
Vm.right.gradient=0;
Vm=setBC(Vm,constField(0,N),xC,xF,0);

% Fluxes initialization
fluxU.internal=zeros(N,1);
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;
fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

% Dummy rho for set BC operations
dummyRho=constField(1,N);
  
% Memory allocation only needed for FVS
if (strcmp(method,'FVS'))
  A=zeros(2,2,N);
end

% Initialization or not
if (initia!=1)
  disp('Continued running');
  load('data.dat');
end

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  printf('Time-step: %d. Time: %g\n',i,i*dt);
  
  % Sum of u along the domain to check conservation
  acc=sum(u.internal);	
  printf('Sum of u in the domain: %g\n',acc);

  % Common fields calculation
  Vpq=assign(constField(V0,N),assign(constField(1,N),u,'-'),'*');
  rhom=assign(assign(constField(rhog,N),u,'*'),assign(assign(constField(1,N),u,'-'),constField(rhol,N),'*'),'+');
  cp=assign(assign(constField(rhog,N),u,'*'),rhom,'/');
  % Vm is calculated from an analytical expression since the complete solver is not available
  Vm=assign(Vpq,assign(u,assign(assign(constField(rhog,N),rhom,'/'),constField(1,N),'-'),'*'),'*');
  % BC's setting
  Vm=setBC(Vm,constField(0,N),xC,xF,0);
	 
  if (strcmp(method,'Rusanov'))

    if (VmIncluded)
      % Vm effect included in total flux
      fluxU=assign(assign(assign(Vpq,u,'*'),assign(constField(1,N),cp,'-'),'*'),assign(Vm,u,'*'),'+');
      % No conservative/centered fluxes are used
      cfluxU=zeros(N+1,1);
    else
      % Vm effect included included as a centered flux in Rusanov scheme
      % In the real solver (rhom*Vm)_f is given by the PISO loop
      if 1
	% Like FOAM
	cfluxU=fvc_interpolate(Vm, w, xC, xF).*fvc_interpolate(u, w, xC, xF);
      else
	cfluxU=fvc_interpolate(assign(Vm,u,'*'),w, xC, xF);
      end
      % Zero centered fluxes at boundaries
      fluxU=assign(assign(Vpq,u,'*'),assign(constField(1,N),cp,'-'),'*');
    end
    
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(u,rhol,rhog,V0,constField(0,N));  
    % Stabilization flux has to be zero at boundaries too
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;  
  
    % Ensure no fluxes at boundaries
    fluxU.left.type='V';
    fluxU.left.value=0;
    fluxU.right.type='V';
    fluxU.right.value=0;
    fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

      
    % Time advancement
    if (VmIncluded)
      [u,u]=Rusanov(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	
    else
      [u]=MarquezNigro(u,u,fluxU,fluxU,cfluxU,cfluxU,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,w,xC,xF,dx,dt);
    end

  elseif (strcmp(method,'Godunov'))

    % No conservative fluxes are used
    cfluxU=zeros(N+1,1);

    [u]=Godunov(u,@alphaEqnVmFlux,cfluxU,dummyRho,dummyRho,V0,rhol,rhog,1,dx,dt,10);

  elseif (strcmp(method,'KT'))

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    % phiAlphag=vanLeer(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9))*0; % Constant values by cells (mimiking Rusanov?)
 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    %[a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,constField(0,N));
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(u,rhol,rhog,V0,0);    
    % Stabilization flux has to be zero at boundaries
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;
    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Time advancement by Kurnanov & Tadmor's scheme
    [uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnVmFluxFlat,@alphaEqnVmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);

  elseif (strcmp(method,'KTcFlux'))

    % Using Kurganov & Tadmor but with Vm as centered flux

    % Limiting
    % Sweby's fuction calculation
    % phiAlphag=superbee(rvalue(u,1E-9));
    phiAlphag=vanLeer(rvalue(u,1E-9));
    %phiAlphag=vanLeer(rvalue(u,1E-9))*0; % Constant values by cells (mimiking Rusanov?)
 
    % Limited values calculation
    [uLimited]=limitedValues(u,phiAlphag,dx,dt);

    if 1  
      % Eigenvalues for complete flux
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(u,rhol,rhog,V0,0);    
    else
      % Eigenvalues for flux without Vm
      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(u,rhol,rhog,V0,0); 
    end
    % Stabilization flux has to be zero at boundaries
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;
    % Flatten for KT function
    aeigens=[a_j_minus_half; a_j_plus_half(end)];

    % Vm as a centered flux
    phic=fvc_interpolate(Vm, w, xC, xF);

    % Impermeable walls
    phic(1)=0;
    phic(end)=0;

    % alphag at faces
    uInt=fvc_interpolate(u, w, xC, xF);

    % Time advancement by Kurnanov & Tadmor's scheme with additional centered flux
    %[uDummy,u]=KT(u,u,uLimited,uLimited,@alphaEqnNoVmFluxFlat,@alphaEqnNoVmFluxFlat,aeigens,dx,dt,V0,rhol,rhog,aexp);
    [u]=KTcFlux(u,uLimited,@alphaEqnNoVmFluxFlat,aeigens,phic,uInt,dx,dt,constField(1,N),constField(1,N),ones(N+1,1),V0,rhol,rhog,aexp);

  else
    % Centered flux for Vm inclusion
    phiAlpha=fvc_interpolate(Vm, w, xC, xF);

    if (strcmp(method,'Centered'))
      % FOAM original
      phiVdrp=fvc_interpolate(assign(Vpq,assign(constField(1,N),cp,'-'),'*'), w, xC, xF);
    elseif (strcmp(method,'UADE'))
    % UADE stabilization
      directionFlux=fvc_interpolate(Vpq, w, xC, xF);
      directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
      phiVdrp=fvc_general_interpolate(assign(Vpq,assign(constField(1,N),cp,'-'),'*'), xC, xF,1,directionFlux);
    end

    phiAlpha=phiVdrp+phiAlpha;

    % Impermeable condition at walls
    phiAlpha(1)=0;
    phiAlpha(end)=0;

    % Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
    directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
    Alphag0Int=fvc_general_interpolate(u, xC, xF,-1,directionFluxAG0);

    % Explicit temporal scheme
    u.internal(1:end)=u.internal(1:end)-dt./dx*(phiAlpha(2:end).*Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1));
  

  end

  % BC adjusting
  u=setBC(u,constField(0,N),xC,xF,g);

  % Printing and saving
  if 1
    hold on;
    if (rem(i,100)==0 || i==1)
      %plot(xC,u.internal,'r*-')
      eval(['save alphaEqnIsolatedVm-' method '-' num2str(i) '.dat u rhom rhom0 Vm Vm0 cp p Vpq Sf xC xF w'])
    end
  end

  % Time debugging variable storing
  TAux(i,1)=sum(u.internal)-acc;
  printf('Delta u in present timestep: %g\n',TAux(i));
  
  % Field actualization
  Vm0=Vm;
  rhom0=rhom;

end

save data.dat u 

