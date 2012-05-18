if 0
  % Methods written in alphag
  if 0
    % UADE-like method
    oneOverRhomF=fvc_interpolate(assign(constField(1,N),rhom0,'/'), w, xC, xF);

    % Phis are velocities
    phiAlpha=rhomPhi.*oneOverRhomF;

    % UADE stabilization
    % Vdrp downwind respect itself
    directionFlux=fvc_interpolate(assign(Vpq,arrayToField(1-cp),'*'), w, xC, xF)+phiAlpha;
    directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
    phiVdrp=fvc_general_interpolate(assign(Vpq,arrayToField(1-cp),'*'), xC, xF,1,directionFlux).*Sf;

    % Final flux
    phiAlpha+=phiVdrp;

    % alphag0 interpolation, upwind respect phiAlpha (gas velocity indeed)
    directionFluxA=sign(phiAlpha(2:end-1,1)+1E-9);
    alphag0Int=fvc_general_interpolate(alphag0, xC, xF,-1,directionFluxA);

    % Impermeable walls
    if 1
      phiAlpha(1)=0;
      %phiAlpha(end)=0;
    end

    % Solve
    if (fullVerbose==1)
      disp('Explicit solving of alphag')
    end	

    % Explicit solving for alphag
    alphag.internal(1:end)=alphag0.internal(1:end)-dt./dx*(phiAlpha(2:end).*alphag0Int(2:end)-phiAlpha(1:end-1).*alphag0Int(1:end-1));
    alphag=setBC(alphag,rhom0,xC,xF,g);

  else

    % Face maximum for spectral radius specially tailored for this problem
    [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnPhiVm(alphag,rhomPhi,rhol,rhog,V0);

    % Stabilization flux has to be zero at boundaries too (impermeable wall)
    a_j_minus_half(1)=0;
    %a_j_plus_half(end)=0;  

    % Left and right face fluxes  
    alphagOverRhom=assign(alphag,rhom0,'/');
    phiVdrp=assign(assign(Vpq,arrayToField(1-cp),'*'),alphag,'*');
    fluxLeft=rhomPhi(1:end-1).*alphagOverRhom.internal(1:end)+phiVdrp.internal(1:end);
    fluxRight=rhomPhi(2:end).*alphagOverRhom.internal(1:end)+phiVdrp.internal(1:end);

    % Centered fluxes
    fluxAlpha=zeros(N+1,1);
    % Interior faces
    fluxAlpha(2:end-1)=(fluxLeft(2:end)+fluxRight(1:end-1))/2;
    % BC's
    fluxAlpha(1)=fluxLeft(1);
    fluxAlpha(end)=fluxRight(end);

    % Impermeable BC's
    if 1
      fluxAlpha(1)=0;
      %fluxAlpha(end)=0;
    end
    
    % Time advancement
    [alphag,alphag]=RusanovB(alphag0,alphag0,fluxAlpha,fluxAlpha,a_j_minus_half,a_j_plus_half,0,0,dx,dt);
    alphag=setBC(alphag,rhom,xC,xF,g);

  end


  % Mixture density prediction by alpha change
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');

  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);

else
  % Methods written in Alphag

  % Local velocities
  [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnPhiVmRhoT(alphag,rhomPhi,rhol,rhog,V0,rhom);

  % Stabilization flux has to be zero at boundaries too (impermeable wall)
  if 1
    a_j_minus_half(1)=0;
    a_j_plus_half(end)=0;  
  end

  % Left and right face fluxes  
  Alphag0=assign(assign(alphag0,constField(rhog,N),'*'),rhom0,'/');
  phiVdrp=assign(assign(Vpq,arrayToField(1-cp),'*'),rhom,'*');
  fluxLeft=(rhomPhi(1:end-1)+phiVdrp.internal).*Alphag0.internal;
  fluxRight=(rhomPhi(2:end)+phiVdrp.internal).*Alphag0.internal;

  % Centered fluxes
  fluxAlpha=zeros(N+1,1);
  % Interior faces
  fluxAlpha(2:end-1)=(fluxLeft(2:end)+fluxRight(1:end-1))/2;
  % BC's
  fluxAlpha(1)=fluxLeft(1);
  fluxAlpha(end)=fluxRight(end);

  % Impermeable BC's
  if 1
    fluxAlpha(1)=0;
    fluxAlpha(end)=0;
  end

  [Alphag,Alphag]=RusanovB(Alphag0,Alphag0,fluxAlpha,fluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,dx,dt);
  Alphag=setBC(Alphag,rhom,xC,xF,g);

  % Mixture density actualization
  rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);

  % alphag actualization
  alphag.internal=rhom.internal.*Alphag.internal./rhog;
  alphag=setBC(alphag,rhom,xC,xF,g);

end


% Save rhom value for debugging purposes
rhomSave=rhom;