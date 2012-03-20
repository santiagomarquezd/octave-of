% Solves alpha advection equation problem by Rosunov's Method
% du/dt+d/dx[F(u)]=0; F(u)=(V0*(1-alphag)*(1-alphag*rhog/rhom)+
% V0*(1-alphag)*alphag*(rhog/rhom-1))*alphag
% The problems is written in terms of A=rhog*alphag/rhom

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
rhog=1;
rhol=1000;
V0=0.282;
alphag0=0.5; %0.5;
g=0; % Dummy for setBC

% Domain extension
xleft=0;
xright=1;

% Cross sectional areas
S=1; %0.02;

% BC's
alphagLeft=0;
alphagRight=1;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=4000;%2000;%0*5000; %100;

% Number of cells
N=100;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;
% Interpolation weights calculation
w=weights(xC, xF);
% Face areas
Sf=ones(size(xF))*S;
% Cell volumes
V=ones(size(xC))*S*dx;

% Fields initialization
% alphag
% alphag=ones(N,1)*alphag0;
alphag.internal=ones(N,1)*alphag0;
alphag.left.type='V';
alphag.left.value=alphagLeft;
alphag.right.type='V';
alphag.right.value=alphagRight;

alphag=setBC(alphag,constField(0,N),xC,xF,0);

% Flux
flux.internal=zeros(N,1);
flux.left.type='V';
flux.left.value=0;
flux.right.type='V';
flux.right.value=0;

dummyRho=constField(1,N);

% rhom field initialization
rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
rhom=setBC(rhom,constField(0,N),xC,xF,0);
rhom0=rhom;

% Alphag initialization;
Alphag.internal=alphag.internal*rhog./rhom.internal;
Alphag.left.type='V';
Alphag.left.value=alphag.left.setvalue*rhog./rhom.left.setvalue;
Alphag.right.type='V';
Alphag.right.value=alphag.right.setvalue*rhog./rhom.right.setvalue;
Alphag=setBC(Alphag,constField(0,N),xC,xF,0);

% Calculation of mean velocity from alphag
Vm.internal=zeros(N,1);%V0.*(1-alphag.internal(1:end)).*alphag.internal(1:end).*(rhog./rhom.internal-1);  
Vm.left.type='V';
Vm.left.value=0;
Vm.right.type='V';
Vm.right.value=0;
Vm=setBC(Vm,constField(0,N),xC,xF,0);


% Calculation of drift velocity from alphag
Vdrp.internal=V0.*(1-alphag.internal(1:end)).*(1-alphag.internal(1:end).*rhog./rhom.internal);  
Vdrp.left.type='V';
Vdrp.left.value=0;
Vdrp.right.type='V';
Vdrp.right.value=0;
Vdrp=setBC(Vdrp,constField(0,N),xC,xF,0);

% Full message selection
fullVerbose=0;  %1: enable, 0:disabled


% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  %i

    % rhomPhi field calculation from Vm and rhom fields
    % Creation of flux direction
    directionFlux=fvc_interpolate(Vm, w, xC, xF);
    directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
    rhomPhi=fvc_interpolate(Vm, w, xC, xF).*fvc_general_interpolate(rhom0, xC, xF,-1,directionFlux).*Sf;

    % Rhom PREDICTOR
    rhoEqn

    %keyboard; pause;

    % Calculation of drift velocity from alphag
    Vdrp.internal=V0.*(1-alphag.internal(1:end)).*(1-alphag.internal(1:end).*rhog./rhom.internal);  
    Vdrp=setBC(Vdrp,constField(0,N),xC,xF,0);  

    %keyboard; pause;
      
    % Rusanov scheme
    % u_j^(n+1)=u_j^n-lambda/2*(f(u_(j+1)^n)-f(u_(j-1)^n))+1/2*(lambda*a_(j+1/2)^n*(u_(j+1)^n-u_j^n)-lambda*a_(j-1/2)^n*(u_j^n-u_(j-1)^n))   

    % Flux evaluation (remember that it is a cell flux, flux at boundaries are stored in BC's)
  
    if 1
      flux=assign(assign(assign(rhom,Vm,'*'),assign(rhom,Vdrp,'*'),'+'),Alphag,'*');
    else
      % Double flux for debugging purposes
      flux=assign(assign(assign(assign(rhom0,Vm,'*'),assign(rhom0,Vdrp,'*'),'+'),Alphag,'*'),constField(2,N),'*');
    end
    
    flux=setBC(flux,constField(0,N),xC,xF,0);

    cflux=zeros(N+1,1);
    
    % Local velocities calculation
    [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag,rhol,rhog,V0);
  
    % Temporal integration by Rusanov scheme
    [Alphag,Alphag]=Rusanov(Alphag,Alphag,flux,flux,cflux,cflux,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,dx,dt);	   

    Alphag=setBC(Alphag,rhom,xC,xF,g);

    % Mixture density CORRECTION
    rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);

    % BC overriden (maybe rhom has to have ZG BC at top and bottom)
    rhom.left.value=rhom.internal(1);
    rhom.right.value=rhom.internal(end);
    rhom=setBC(rhom,constField(0,N),xC,xF,g);

    % alpha actualization
    alphag.internal=rhom.internal.*Alphag.internal./rhog;
    alphag=setBC(alphag,rhom,xC,xF,g);
  
    if 0
  
      figure(1); plot(xC,alphag);
      figure(2); plot(eig1);
      figure(3); plot(eig2);
      figure(4); plot(xC,Vm);
  
    end

    % Velocities actualization
    Vm.internal=V0.*(1-alphag.internal(1:end)).*alphag.internal(1:end).*(rhog./rhom0.internal-1);  
    Vm=setBC(Vm,constField(0,N),xC,xF,0);

    rhom0=rhom;

end
