% -----------------------------------------
% Explicit solution of Alpha by UADE method

if 0
  % This portion of code is not present in FOAM original

  % Creates Alphag0 and Alphag as a copy of alpha.
  Alphag0=alphag;
  Alphag=Alphag0;

  % Determination of Alpha auxiliary field
  %Alphag0=assign(assign(alphag0,constField(rhog,N),'*'),rhom,'/');
  Alphag0.internal=alphag0.internal.*rhog./rhom0.internal;
  Alphag0=setBC(Alphag0,rhom0,xC,xF,g);
end

% alphaEqn
phiAlpha=rhomPhi;

% Final flux, PISO part
fluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF);

% UADE stabilization with full upwind/downwind
directionFlux=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF)+phiAlpha;
directionFlux=sign(directionFlux(2:end-1,1)+1E-9);
phiVdrp=fvc_general_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), xC, xF,1,directionFlux).*Sf;

phiAlpha+=phiVdrp;

% Alphag0 values at interfaces with full upwind (direction given by phiAlpha)
directionFluxAG0=sign(phiAlpha(2:end-1,1)+1E-9);
Alphag0Int=fvc_general_interpolate(Alphag0, xC, xF,-1,directionFluxAG0);

% Final flux, relative velocity part
fluxAlpha=fluxAlpha+phiVdrp.*Alphag0Int;

% Solve
if (fullVerbose==1)
  disp('Explicit solving of Alphag')
end	

if 0
  % Explicit solving for Alphag
  Alphag.internal(1:end)=Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx*(phiAlpha(2:end).*...
			  Alphag0Int(2:end)-phiAlpha(1:end-1).*Alphag0Int(1:end-1))./rhom.internal(1:end);
else
  % Explicit solving for Alphag with consolided flux
  
  if 0
	% Impermeable walls test
	phiAlpha(1)=0;
	phiAlpha(end)=0;
  end
  
  % Explicit solving for AlphagRhom
  if 1
    % Overwrites the flux with a new one where Alpha is ever Downwind/Upwind
    fluxAlpha=phiAlpha.*Alphag0Int;
  end

  Alphag.internal(1:end)=Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx*(fluxAlpha(2:end)-...
			  fluxAlpha(1:end-1))./rhom.internal(1:end);
end

% -----------------------------------------
% Explicit solution of Alpha by Rusanov method

if 0
%    %Old version (working)
%  
%    % Face maximum for spectral radius specially tailored for this problem
%    if 0
%      % U-alphag block
%      [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag0,rhol,rhog,V0);
%    else
%      % Isolated alpha equation with explicit relation for Vm (not calculating dF/dalphag)
%      [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(alphag,rhol,rhog,V0);
%    end
%  
%    if 1
%      % Stabilization flux has to be zero at boundaries too (impermeable wall)
%      a_j_minus_half(1)=0;
%      %a_j_plus_half(end)=0;  
%    end
%  
%    % Quasi-conservative, rhomPhi is used instead of creating it
%    fluxVdrpRF=assign(assign(rhom0,assign(Vpq,arrayToField(1-cp),'*'),'*'),Alphag0,'*');
%    fluxVdrpRF=setBC(fluxVdrpRF,constField(0,N),xC,xF,0);
%  
%    if 0
%      %Ensure no fluxes at boundaries
%      fluxVdrpRF.left.type='V';
%      fluxVdrpRF.left.value=0;
%      fluxVdrpRF.right.type='V';
%      fluxVdrpRF.right.value=0;
%      fluxVdrpRF=setBC(fluxAlpha,constField(0,N),xC,xF,0);
%    end
%  
%    if 0
%      cFluxAlpha=rhomPhi.*fvc_interpolate(Alphag0, w, xC, xF);
%      % cFluxAlpha=fvc_interpolate(assign(assign(rhom,Vm,'*'),Alphag0,'*'), w, xC, xF);	
%      [AlphagRF,AlphagRF]=Rusanov(Alphag0,Alphag0,fluxVdrpRF,fluxVdrpRF,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,dx,dt);  
%    elseif 0
%      cFluxAlpha=rhomPhi;
%      [AlphagRF]=MarquezNigro(Alphag0,Alphag0,fluxVdrpRF,fluxVdrpRF,cFluxAlpha,cFluxAlpha,a_j_minus_half,a_j_plus_half,0,0,rhom0,rhom,w,xC,xF,dx,dt);
%    else
%      % Manual solution taking values at faces
%      % Flux calculation
%      % Internal faces
%      fluxVdrpRFf=zeros(N+1,1);
%      % Standard central difference  
%      fluxVdrpRFf(2:end-1)=(fluxVdrpRF.internal(1:end-1)+fluxVdrpRF.internal(2:end))/2;
%  
%      % Stabilization
%      % It requires the presence of rhom0, so it is interpolated at faces
%      if 1
%        rhom0f=fvc_interpolate(rhom0, w, xC, xF);
%      elseif 0
%        rhom0f=fvc_interpolate(rhom, w, xC, xF);
%      else
%        % Direction given by fluxVdrpRFf, do upwind
%        directionFlux=sign(fluxVdrpRFf(2:end-1,1)+1E-9);
%        rhom0f=fvc_general_interpolate(rhom0, xC, xF,1,directionFlux);
%      end
%      % Final flux, internal faces
%      % Using eigenvalues based in alpha equation (not A) leads to rhom0 multiplication
%      % Here it'd be used eigenvalues calculated as dF/dA like in aspeesIsolatedAlphaEqnPhiVmRhoT.m
%      % it has into account the frozen values of Vm*rhom given in rhomPhi and the eventual frozen value of 
%      % rhom prediction rhomT. Both of these values, rhomPhi and rhonT (rhom) are given.  
%      fluxVdrpRFf(2:end-1)=fluxVdrpRFf(2:end-1)-a_j_plus_half(1:end-1).*rhom0f(2:end-1).*(Alphag0.internal(2:end)-Alphag0.internal(1:end-1))/2;
%      % Boundary faces
%      fluxVdrpRFf(1)=fluxVdrpRF.left.setvalue+a_j_minus_half(1).*(Alphag0.internal(1)-Alphag0.left.setvalue)*rhom0f(1);
%      fluxVdrpRFf(end)=fluxVdrpRF.right.setvalue+a_j_plus_half(1).*(Alphag0.right.setvalue-Alphag0.internal(end))*rhom0f(end);
%      % Flux from PISO is added multiplied by the corresponding Alphag0 values
%      Alphag0f=fvc_interpolate(Alphag0, w, xC, xF);
%  
%      % Semi-analytic version of flux
%      % Vm is calculated from an analytical expression
%      Vm=assign(Vpq,assign(alphag,assign(assign(constField(rhog,N),rhom0,'/'),constField(1,N),'-'),'*'),'*');
%      fluxVdrpRFAnalytic=assign(assign(rhom0,assign(Vpq,arrayToField(1-cp),'*'),'*'),Alphag0,'*');
%      fluxAlphaRFAnalytic=assign(assign(assign(rhom0,Vm,'*'),Alphag0,'*'),fluxVdrpRFAnalytic,'+');
%  
%      if 1
%        % Use central difference for Alphag0f
%        fluxAlphaRFf=fluxVdrpRFf+rhomPhi.*Alphag0f;
%      elseif 0
%        % Use upwinding for Alphag0f, this value is given by Alphag0Int calculated in UADE's section
%        fluxAlphaRFf=fluxVdrpRFf+rhomPhi.*Alphag0Int;
%      elseif 0
%        % Flux for mixture velocity is not taken from PISO
%        fluxAlphaRFf=fluxVdrpRFf+fvc_interpolate(assign(assign(rhom,U,'*'),Alphag0,'*'), w, xC, xF);
%      else
%        % Use semi-analytic flux
%        fluxAlphaRFf=fvc_interpolate(fluxAlphaRFAnalytic, w, xC, xF);
%      end
%      
%  
%      % Time advancement
%      % Memory allocation
%      AlphagRF=Alphag0;
%      % Calculation
%      AlphagRF.internal(1:end)=Alphag0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-dt./dx./rhom.internal(1:end).*(fluxAlphaRFf(2:end)-fluxAlphaRFf(1:end-1));
%      AlphagRF=setBC(AlphagRF,rhom,xC,xF,g);
%    end

else

  % New version based in K&T scheme, High Resolution schemes can be activated (testing)
  % ****************************************************************************  
  % IMPORTANT: Flux calculation and stabilization is based directly on alphag
  % ****************************************************************************

  % Limiting
  % Sweby's fuction calculation
  % phiAlphag=superbee(rvalue(alphag0,1E-9));
  % phiAlphag=vanLeer(rvalue(alphag0,1E-9));
   phiAlphag=vanLeer(rvalue(alphag0,1E-9))*0; % Constant values by cells (mimiking Rusanov?)

  % Limited values calculation
  [alphagLimited]=limitedValues(alphag0,phiAlphag,dx,dt);

  % Here it'd be used eigenvalues calculated as dF/dA like in aspeesIsolatedAlphaEqnPhiVmRhoT.m
  % it has into account the frozen values of Vm*rhom given in rhomPhi and the eventual frozen value of 
  % rhom prediction rhomT. Both of these values, rhomPhi and rhonT (rhom) are given.  
  [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVm(alphag0,rhol,rhog,V0,aexp);    
  %[a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnVmGastaldo(alphag,rhol,rhog,V0); % Gastaldo's example
  % Stabilization flux has to be zero at boundaries
  a_j_minus_half(1)=0;
  a_j_plus_half(end)=0;
  % Flatten for KT function
  aeigens=[a_j_minus_half; a_j_plus_half(end)];

  % Alpha at faces
  AlphagInt=fvc_interpolate(Alphag0, w, xC, xF);

  % rhom0 at faces
  rhom0Int=fvc_interpolate(rhom0, w, xC, xF);

  % Time advancement by Kurnanov & Tadmor's scheme
  [AlphagRF]=KTcFlux(Alphag0,alphagLimited,@AEqnNoVmFluxFlat,aeigens,rhomPhi,AlphagInt,dx,dt,rhom,rhom0,rhom0Int,V0,rhol,rhog,aexp);

  %keyboard; pause;
  %disp('Hola')

end

% -----------------------------------------
% Temporal actualization of remaining variables (rhom and alphag)

disp('A antes de actualizar')
sum(Alphag0.internal)
disp('A despues de actualizar (UADE)')
sum(Alphag.internal)
disp('A despues de actualizar (Rusanov)')
sum(AlphagRF.internal)

% Uses UADE solution (0) or Rusanov (1)
if 1
  Alphag=AlphagRF;
end

% BC setting
Alphag=setBC(Alphag,rhom,xC,xF,g);

% Mixture density actualization
rhom.internal=rhol./(1+(rhol/rhog - 1.0).*Alphag.internal);
% BC overriden (maybe rhom has to have ZG BC at top and bottom)
rhom.left.value=rhom.internal(1);
rhom.right.value=rhom.internal(end);
rhom=setBC(rhom,constField(0,N),xC,xF,g);

% rhom adjustement in first step
if (step==1)
  factor=sum(rhom0.internal)/sum(rhom.internal)
  rhom.internal=rhom.internal*factor;
  % BC overriden (maybe rhom has to have ZG BC at top and bottom)
  rhom.left.value=rhom.internal(1);
  rhom.right.value=rhom.internal(end);
  rhom=setBC(rhom,constField(0,N),xC,xF,g);
end

% alpha actualization
alphag.internal=rhom.internal.*Alphag.internal./rhog;
alphag=setBC(alphag,rhom,xC,xF,g);

rhomSave=rhom;

% UADE
%  flux1=phiVdrp(1:end-1).*Alphag0Int(1:end-1);
%  flux2=fluxAlphaRFf;
%  % flux2 with stabilization
%  flux3=(flux2.internal(3:end)+flux2.internal(2:end-1))/2+a_j_plus_half(2:end-1).*(Alphag0.internal(3:end)-Alphag0.internal(2:end-1));
%  %
%  flux4=fvc_interpolate(assign(rhom,assign(Vpq,arrayToField(1-cp),'*'),'*'), w, xC, xF).*fvc_interpolate(Alphag, w, xC, xF);
%  flux4=flux4(2:end-1)+a_j_plus_half(1:end-1).*(Alphag0.internal(2:end)-Alphag0.internal(1:end-1));
%  

