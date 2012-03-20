  % Face maximum for spectral radius specially tailored for this problem
  [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag0,rhol,rhog,V0);

  % alphaFlux
  % Calculation of drift velocity from alphag
  Vdrp=U;
  Vdrp.internal=Vpq.internal(1:end).*(1-alphag0.internal(1:end).*rhog./rhom.internal);  
  Vdrp=setBC(Vdrp,constField(0,N),xC,xF,0);  

  % Flux evaluation (remember that it is a cell flux, flux at boundaries are stored in BC's)
  fluxAlpha=assign(assign(U0,Vdrp,'+'),alphag0,'*');
  fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);
  cFluxAlpha=zeros(N+1,1);

  % U flux
  fluxU=assign(U0,assign(U0,rhom,'*'),'*');
  fluxU=setBC(fluxU,constField(0,N),xC,xF,0);
  cFluxU=zeros(N+1,1);

  % U source term
  % Drift term
  arg=assign(assign(rhom,arrayToField((1-cp).*cp),'*'),assign(Vpq,constField(2,N),'^'),'*');
  drift=-fvc_div_cell(arg, w, xC, xF, Sf, V);

  if 1
  
    % UEqn implicit assembling for PISO solver
    [ddtM, ddtRHS] = fvm_ddt(rhom,rhom0,U0,V,dt,1);
    [convM, convRHS]= fvm_div_flux_cell(rhomPhi,U0,xC,xF,w,1);

    % Calculating explicit terms (to RHS)
    driftRHS=drift.*V;

    % Final assembling
    UEqnM=ddtM+convM;
    UEqnRHS=ddtRHS+convRHS+driftRHS;

  end


  % Gravitational and pressure terms
  vol=-gh.*fvc_grad(rhom, w, xC, xF, Sf, V)-fvc_grad(p_rgh, w, xC, xF, Sf, V);
  SU=drift+vol;

  if 1
    % Non-conservative
    [U,alphag]=Rusanov(U0,alphag0,fluxU,fluxAlpha,cFluxU,cFluxAlpha,a_j_minus_half,a_j_plus_half,SU,0,rhom0,rhom,dx,dt);
    
    % U boundary setting (the obtained velocity is the momentum predictor of the problem)
    U=setBC(U,rhom,xC,xF,g);
    % alphag boundary setting
    alphag=setBC(alphag,rhom,xC,xF,g);
  elseif 1
    %
  else
    %
  end

  % Secondary phase bounding
  if 0
    alphag=bound(alphag,'min',0);
    alphag=bound(alphag,'max',1);
    alphag=setBC(alphag,rhom,xC,xF,g);
  end

  % rhom field actualization
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
  rhom=setBC(rhom,constField(0,N),xC,xF,0);

  % For explicit PISO solver
  Ap=rhom.internal/dt;
  