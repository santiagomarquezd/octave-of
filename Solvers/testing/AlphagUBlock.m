  % The advected quantity for alpha equation is Alpha  
  % Alllocation
  Alphag0=alphag;
  % Value setting
  Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
  Alphag0=setBC(Alphag0,rhom,xC,xF,g);   

  % Face maximum for spectral radius specially tailored for this problem
  [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag0,rhol,rhog,V0);

  % AlphaFlux
  fluxAlpha=assign(assign(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'),assign(U,rhom,'*'),'+'),Alphag0,'*');
  fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);

  % U flux
  fluxU=assign(U0,assign(U0,rhom,'*'),'*');
  fluxU=setBC(fluxU,constField(0,N),xC,xF,0);

  % U source term
  % Drift term
  arg=assign(assign(rhom,arrayToField((1-cp).*cp),'*'),assign(Vpq,constField(2,N),'^'),'*');
  drift=-fvc_div_cell(arg, w, xC, xF, Sf, V);

  % Gravitational and pressure terms
  vol=-gh.*fvc_grad(rhom, w, xC, xF, Sf, V)-fvc_grad(p_rgh, w, xC, xF, Sf, V);
  SU=drift+vol;

  %Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx*(phiAlpha(2:end).*...
  % 			 Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1))./rhom.internal(1:end);

  % Time advancement
  if 1
    % Non-conservative
    cFluxAlpha=rhomPhi*0;
    cFluxU=rhomPhi*0;
    [U,Alphag]=Rusanov(U0,Alphag0,fluxU,fluxAlpha,cFluxU,cFluxAlpha,a_j_minus_half,a_j_plus_half,SU,0,rhom0,rhom,dx,dt);
  elseif 1
    % Quasi-conservative, rhomPhi is used instead of creating it

    % AlphaFlux (only relative velocity part, Um flux is taken from PISO)
    if 0
      % Overrides rhomPhi flux from PISO by a local version based on Vm(alphag) calculation
      % Full centered flux
      Vm=U;
      Vm.internal=V0.*(1-alphag.internal(1:end)).*alphag.internal(1:end).*(rhog./rhom.internal-1);
      Vm=setBC(Vm,constField(0,N),xC,xF,0);
      %rhomPhi=fvc_interpolate(assign(rhom,Vm,'*'), w, xC, xF);	
      %rhomPhi=fvc_interpolate(rhom, w, xC, xF).*fvc_interpolate(Vm, w, xC, xF);	
    end

    fluxAlpha=assign(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'),Alphag0,'*');
    fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);
    
    %cFluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF)*0;
    %cFluxAlpha=fvc_interpolate(assign(assign(rhom,Vm,'*'),Alphag0,'*'), w, xC, xF);	
    
    % U flux [(rhom*U)_f*Sf_*interpolate(U)]. All the flux is in cFluxU
    % Zeroing fluxU
    fluxU=assign(fluxU,constField(0,N),'*');
    cFluxU=rhomPhi.*fvc_interpolate(U0, w, xC, xF);

    [U,Alphag]=Rusanov(U0,Alphag0,fluxU,fluxAlpha,cFluxU,cFluxAlpha,a_j_minus_half,a_j_plus_half,SU,0,rhom0,rhom,dx,dt);  

  else
    % Conservative for U and nonconservative for alpha
    cFluxAlpha=rhomPhi*0;  
    fluxU=assign(fluxU,constField(0,N),'*');
    cFluxU=rhomPhi.*fvc_interpolate(U0, w, xC, xF);

    [U,Alphag]=Rusanov(U0,Alphag0,fluxU,fluxAlpha,cFluxU,cFluxAlpha,a_j_minus_half,a_j_plus_half,SU,0,rhom0,rhom,dx,dt);  

  end
	  

  if 1
  Alphag=setBC(Alphag,rhom,xC,xF,g);

  % Mixture density actualization
  rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);

  % alpha actualization
  alphag.internal=rhom.internal.*Alphag.internal./rhog;
  alphag=setBC(alphag,rhom,xC,xF,g);
  else
    alphaEqnIshii
  end



  % Secondary phase bounding
  if 0
    alphag=bound(alphag,'min',0);
    alphag=bound(alphag,'max',1);
    alphag=setBC(alphag,rhom,xC,xF,g);
  end

  %keyboard; pause; 

  if 1

  % U boundary setting (the obtained velocity is the momentum predictor of the problem)
  U=setBC(U,rhom,xC,xF,g);

  % UEqn implicit assembling for PISO solver
  [ddtM, ddtRHS] = fvm_ddt(rhom,rhom0,U0,V,dt,1);
  [convM, convRHS]= fvm_div_flux_cell(rhomPhi,U0,xC,xF,w,1);

  % Calculating explicit terms (to RHS)
  driftRHS=drift.*V;

  % Final assembling
  UEqnM=ddtM+convM;
  UEqnRHS=ddtRHS+convRHS+driftRHS;

  end

  % For explicit PISO solver
  Ap=rhom.internal/dt;
  