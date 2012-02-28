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
  %phiVdrp=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF).*Sf;
  directionFlux=fvc_interpolate(Vpq, w, xC, xF);
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;
  phiAlpha+=phiVdrp;
else
  directionFlux=fvc_interpolate(Vg, w, xC, xF);
  directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
  phiVdrp=fvc_general_interpolate(assign(rhom,Vg,'*'), xC, xF,1,directionFlux).*Sf;
  phiAlpha+=phiVdrp;
end


% Assembling
if (fullVerbose==1)
  disp('Assembling AlphaGEqn')
end
[ddtMat, ddtRHS] = fvm_ddt(rhom,rhom0,Alphag0,V,dt,1);
[divMat, divRHS]= fvm_div_flux_cell(phiAlpha,Alphag0,xC,xF,w,alphaDiv);
[lapMat, lapRHS]= fvm_laplacian(Alphag0,nu,xC,xF,Sf);
  
%  M=M+A+B;
%  RHS=RHS+ARHS+BRHS;
AlphaM=ddtMat+divMat-lapMat;
RHS=ddtRHS+divRHS-lapRHS;
 
% Solve
if (fullVerbose==1)
  disp('Solving for Alphag')
end
Alphag.internal=AlphaM\RHS;
Alphag=setBC(Alphag,rhom,xC,xF,g);



%  flux.internal=phiAlpha(2:end-1);
%  flux.left.setvalue=phiAlpha(1);
%  flux.right.setvalue=phiAlpha(end);

%keyboard; pause;

% Alphag.internal=LxF(Alphag0,Alphag0,flux,flux,xC(2)-xC(1),dt)


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

% rhomPhi field calculation from U and rhom fields
% Creation of flux direction
directionFlux=fvc_interpolate(U, w, xC, xF);
directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
rhomPhi=fvc_interpolate(U, w, xC, xF).*fvc_general_interpolate(rhom, xC, xF,-1,directionFlux).*Sf;
%  rhomPhi=fvc_interpolate(U, w, xC, xF).*fvc_interpolate(rhom, w, xC, xF).*Sf;
