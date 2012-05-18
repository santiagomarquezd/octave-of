  % The advected quantity for alpha equation is Alpha  
  % Alllocation
  Alphag0=alphag;
  % Value setting
  Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
  Alphag0=setBC(Alphag0,rhom,xC,xF,g);   

  % Matrix for the quasi-linear system of U-Alpha
  for i=1:N
    Aarray(:,:,i)=[2*U0.internal(i)-(alphag.internal(i)*(rhog-rhol)*U0.internal(i))/rhom0.internal(i),...
    ((rhog-rhol)*U0.internal(i)^2)/rhom0.internal(i)-((rhog-rhol)*U0.internal(i)*...
    (alphag.internal(i)*(-(2*(1-alphag.internal(i))*rhol*V0)/((1-alphag.internal(i))*...
    rhol+alphag.internal(i)*rhog)-((1-alphag.internal(i))^2*(rhog-rhol)*rhol*V0)/...
    ((1-alphag.internal(i))*rhol+alphag.internal(i)*rhog)^2)+((1-alphag.internal(i))^2*rhol*V0)/...
    ((1-alphag.internal(i))*rhol+alphag.internal(i)*rhog)+U0.internal(i)))/rhom0.internal(i);
    alphag.internal(i),alphag.internal(i)*(-(2*(1-alphag.internal(i))*rhol*V0)/...
    ((1-alphag.internal(i))*rhol+alphag.internal(i)*rhog)-((1-alphag.internal(i))^2*...
    (rhog-rhol)*rhol*V0)/((1-alphag.internal(i))*rhol+alphag.internal(i)*rhog)^2)+...
    ((1-alphag.internal(i))^2*rhol*V0)/((1-alphag.internal(i))*rhol+alphag.internal(i)*rhog)+U0.internal(i)];
  end
  
  % Eigenvalues and eigenvectors for A matrices
  [v,vT,lambda]=arrayEig(Aarray);

  % A matrix spllittings
  % Aplus=matrixArrayProd(matrixArrayProd(v,lambdaPlus(lambda)),vT);
  % Aminus=matrixArrayProd(matrixArrayProd(v,lambdaMinus(lambda)),vT);
  
  % Absolute value of A
  [AAbs]=matrixArrayAbs(v,lambda);  

  % Fluxes calculation
  % Part from PISO (by face)
  cFluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF);
  
  % Part from Drift (flux by cell)
  fluxAlpha=assign(assign(assign(rhom0,assign(Vpq,arrayToField(1-cp),'*'),'*'),assign(U0,rhom0,'*'),'+'),Alphag0,'*');
  
  % Time advancement
  [Alphag]=FVScFlux(U0,Alphag0,AAbs,fluxAlpha,cFluxAlpha,rhom,rhom0,dx,dt,N,xC,xF,w);

  % Apply BC	    
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
  