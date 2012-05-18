% Face maximum for spectral radius specially tailored for this problem
if 0
  % U-alphag block
  [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag0,rhol,rhog,V0);
else
  % Isolated alpha equaation
  [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(alphag0,rhol,rhog,V0,U);
end

% Stabilization flux has to be zero at boundaries too (impermeable wall)
a_j_minus_half(1)=0;
a_j_plus_half(end)=0;  

% Time advancement
% Fully-conservative non using the flux given by PISO
% alphaFlux
if 0
  % Advecting with predicted U (by Chorin's method)
  fluxAlpha=assign(assign(assign(Vpq,alphag0,'*'),assign(constField(1,N),arrayToField(cp),'-'),'*'),assign(U,alphag0,'*'),'+');
else
  % adventing with U0, consistent with alphag0
  fluxAlpha=assign(assign(assign(Vpq,alphag0,'*'),assign(constField(1,N),arrayToField(cp),'-'),'*'),assign(U0,alphag0,'*'),'+');
end

% Vg calculation for debugging purposes
Vg=assign(assign(Vpq,assign(constField(1,N),arrayToField(cp),'-'),'*'),U0,'+');

% Ensure no fluxes at boundaries
fluxU.left.type='V';
fluxU.left.value=0;
fluxU.right.type='V';
fluxU.right.value=0;
fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

cFluxAlpha=rhomPhi*0;
% Dummy rhoField for function compatibility
dummyRho=constField(1,N);

% Time advancement
[alphag,alphag]=Rusanov(alphag0,alphag0,fluxAlpha,fluxAlpha,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);
alphag=setBC(alphag,rhom,xC,xF,g);

% Mixture density prediction by alpha change
rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');

% BC overriden (maybe rhom has to have ZG BC at top and bottom)
rhom.left.value=rhom.internal(1);
rhom.right.value=rhom.internal(end);
rhom=setBC(rhom,constField(0,N),xC,xF,g);


% Save rhom value for debuggin purposes
rhomSave=rhom;
