% Pressure post-processing from Vm, rhom
% Model parameters
g=-10;

% Method selection, Rusanov, Godunov, Centered, UADE
method='Rusanov'; %'KTcFlux'; %'Rusanov';

% Cell number
N=400; %1000;

% End step
nSteps=1000;
step=100;

% Dynamic p
pDynF=zeros(N+1,1);

% File sequence generation
seq=[1 step:step:nSteps];

for step=seq
  % Prints the actual iteration
  printf('Time-step: %d.\n',step);

  % Loads file
  eval(['load alphaEqnIsolatedVm-' method '-' num2str(step) '.dat']);

  % Cell volumes
  V=ones(size(xC))*Sf(1)*dx;

  % Pressure calculation via momentum equation intregation
  % Volumetric term
  temporal=fvc_ddt(assign(Vm,rhom,'*'), assign(Vm0,rhom0,'*'), dt);
  arg=assign(assign(rhom,assign(assign(constField(1,N),cp,'-'),cp,'*'),'*'),assign(Vpq,constField(2,N),'^'),'*');
  drift=fvc_div_cell(arg, w, xC, xF, Sf, V);
  gravity=assign(constField(-g,N),rhom,'*');
  volumetric=assign(assign(assign(arrayToField(temporal),arrayToField(drift),'+'),gravity,'+'),arrayToField(V),'*');

  % Value at top (BC)
  pDynF(N+1)=rhom.right.setvalue*Vm.right.setvalue^2;

  % Pressure integration at faces from top values
  for i=N:-1:1
    pDynF(i)=pDynF(i+1)+volumetric.internal(i); 
  end

  % Static p extraction
  % pF=pDynF-(rhom*Vm^2)_f
  % Momentum at faces estimation
  rhomVm2F=fvc_interpolate(assign(rhom,assign(Vm,Vm,'*'),'*'), w, xC, xF);
  pF=pDynF-rhomVm2F;

  % Pressure field allocation
  pS.internal=ones(N,1)*0;
  pS.left.type='BP';
  pS.right.type='V';
  pS.right.value=0;
  pS=setBC(pS,rhom,xC,xF,g);

  % Pressure reconstruction at cell centers
  pS.internal=fvc_reconstruct(pF,Sf);
  pS=setBC(pS,rhom,xC,xF,g);

  % Saves pressure file
  eval(['save alphaEqnIsolatedVmP-' method '-' num2str(step) '.dat xF xC pF pS'])

end

