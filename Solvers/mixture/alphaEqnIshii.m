% Numerical diffusivity
nu=mult*1/2*mean(abs(U.internal)+abs(Vpq.internal))*mean(dx)*ones(size(rhomPhi));

% nu overriden for test purposes
if 1
  nu=0.00705*ones(size(rhomPhi)); %00705
end

% Creates Alphag0 and Alphag as a copy of alpha.
Alphag0=alphag;
Alphag=Alphag0;

% Determination of Alpha auxiliary field
%Alphag0=assign(assign(alphag0,constField(rhog,N),'*'),rhom,'/');
Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
Alphag0=setBC(Alphag0,rhom,xC,xF,g);

% alphaEqn
phiAlpha=rhomPhi;
if 1
  % UADE stabilization
  directionFlux=fvc_interpolate(Vpq, w, xC, xF);
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;
else 
  % UADE stabilization with full upwind/downwind
  directionFlux=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF)+phiAlpha;
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;
end
phiAlpha+=phiVdrp;

if 0
  % Impermeable walls test
  phiAlpha(1)=0;
  phiAlpha(end)=0;
end


if (alphaExplicit==0)

  % Assembling
  if (fullVerbose==1)
    disp('Assembling AlphaGEqn')
  end
  [ddtMat, ddtRHS] = fvm_ddt(rhom,rhom0,Alphag0,V,dt,1);
  [divMat, divRHS]= fvm_div_flux_cell(phiAlpha,Alphag0,xC,xF,w,alphaDiv);
  [lapMat, lapRHS]= fvm_laplacian(Alphag0,nu,xC,xF,Sf);

  % Solve
  if (fullVerbose==1)
    disp('Implicit solving of Alphag')
  end
  %  M=M+A+B;
  %  RHS=RHS+ARHS+BRHS;
  AlphaM=ddtMat+divMat-lapMat;
  RHS=ddtRHS+divRHS-lapRHS;
  Alphag.internal=AlphaM\RHS;

  %keyboard; pause;

elseif (alphaExplicit==1)
  
  % Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
  directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
  Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionFluxAG0);

  % Solve
  if (fullVerbose==1)
    disp('Explicit solving of Alphag')
  end	

  % Explicit solving for Alphag
  Alphag.internal(1:end)=Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx*(phiAlpha(2:end).*...
  			 Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1))./rhom.internal(1:end);


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

%end
else 
  alphag=setBC(alphag,rhom,xC,xF,g);

  % Mixture density actualization
  rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);
end

%load densities.dat
rhomSave=rhom;