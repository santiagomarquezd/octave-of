  % The advected quantity for alpha equation is Alpha  
  % Alllocation
  Alphag0=alphag;
  % Value setting
  Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
  Alphag0=setBC(Alphag0,rhom,xC,xF,g);   

  % Face maximum for spectral radius specially tailored for this problem
  if 1
    % U-alphag block
    [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag0,rhol,rhog,V0);
  else
    % Isolated alpha equaation
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqn(alphag,rhol,rhog,V0,U);
  end

  % Stabilization flux has to be zero at boundaries too (impermeable wall)
  a_j_minus_half(1)=0;
  a_j_plus_half(end)=0;  

  % Time advancement
  if 0
    % Fully-conservative non using the flux given by PISO
    % AlphaFlux
    fluxAlpha=assign(assign(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'),assign(U,rhom,'*'),'+'),Alphag0,'*');
    fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);
    cFluxAlpha=rhomPhi*0;
    [Alphag,Alphag]=Rusanov(Alphag0,Alphag0,fluxAlpha,fluxAlpha,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,dx,dt);
  elseif 1
    % Quasi-conservative, rhomPhi is used instead of creating it
    fluxAlpha=assign(assign(rhom0,assign(Vpq,arrayToField(1-cp),'*'),'*'),Alphag0,'*');
    fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);

    % Ensure no fluxes at boundaries
%      fluxAlpha.left.type='V';
%      fluxAlpha.left.value=0;
%      fluxAlpha.right.type='V';
%      fluxAlpha.right.value=0;
%      fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);

    if 0
      cFluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF);
      % cFluxAlpha=fvc_interpolate(assign(assign(rhom,Vm,'*'),Alphag0,'*'), w, xC, xF);	
      [Alphag,Alphag]=Rusanov(Alphag0,Alphag0,fluxAlpha,fluxAlpha,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,dx,dt);  
    else
      cFluxAlpha=rhomPhi;
      [Alphag]=MarquezNigro(Alphag0,Alphag0,fluxAlpha,fluxAlpha,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,w,xC,xF,dx,dt);
    end

  else
    %disp('Hola')
    % Non conservative downwind flux for Vdrp (the main idea is testing wheter the problem is in U stabilization or Vdrp)
    fluxAlpha=assign(assign(U,rhom,'*'),Alphag0,'*');
    fluxAlpha=setBC(fluxAlpha,constField(0,N),xC,xF,0);
    directionFlux=fvc_interpolate(Vpq, w, xC, xF);
    directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
    cFluxAlpha=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;
    [Alphag,Alphag]=Rusanov(Alphag0,Alphag0,fluxAlpha,fluxAlpha,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,dx,dt);
  end
	  
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

  % Secondary phase bounding
  if 0
    alphag=bound(alphag,'min',0);
    alphag=bound(alphag,'max',1);
    alphag=setBC(alphag,rhom,xC,xF,g);
  end

  %load density.dat
  rhomSave=rhom;
  