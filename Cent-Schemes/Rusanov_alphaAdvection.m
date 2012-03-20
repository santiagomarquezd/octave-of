% Solves alpha advection equation problem by Rosunov's Method
% du/dt+d/dx[F(u)]=0; F(u)=(V0*(1-alphag)*(1-alphag*rhog/rhom)+
% V0*(1-alphag)*alphag*(rhog/rhom-1))*alphag

% Variables clearance
clear all;
%close all;
page_screen_output(0);

% Physical paramaters
rhog=1;
rhol=1000;
V0=0.282;
alphag0=0.5; %0.5;

% Domain extension
xleft=0;
xright=1;


% BC's
alphagLeft=0;
alphagRight=1;

% Time-step
dt=0.001;

% Number of timesteps
timesteps=1; %100;

% Number of cells
N=100;

% Numerical Pre-processing

% Cell centers and face centers calculation
% Equi-spaced cells
dx=(xright-xleft)/N;
xC=((xleft+dx/2):dx:(xright-dx/2))';
xF=(xleft:dx:xright)';
lambda=dt/dx;

% Fields initialization
% alphag
% alphag=ones(N,1)*alphag0;
alphag.internal=ones(N,1)*alphag0;
alphag.left.type='V';
alphag.left.value=alphagLeft;
alphag.right.type='V';
alphag.right.value=alphagRight;

alphag=setBC(alphag,constField(0,N),xC,xF,0);

% Flux
flux.internal=zeros(N,1);
flux.left.type='V';
flux.left.value=0;
flux.right.type='V';
flux.right.value=0;

dummyRho=constField(1,N);

% Vm and Vdrp allocation
Vm=flux;
Vdrp=flux;

% Temporal loop
for i=1:timesteps

  % Prints the actual iteration
  %i
  
    if 0
      % Rusanov scheme
      % u_j^(n+1)=u_j^n-lambda/2*(f(u_(j+1)^n)-f(u_(j-1)^n))+1/2*(lambda*a_(j+1/2)^n*(u_(j+1)^n-u_j^n)-lambda*a_(j-1/2)^n*(u_j^n-u_(j-1)^n))   

      % Flux evaluation (remember that it is a cell flux, flux at boundaries are stored in BC's)
      
      % rhom field
      rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
      rhom=setBC(rhom,constField(0,N),xC,xF,0);

      flux.internal=(V0.*(1-alphag.internal(1:end)).*(1-alphag.internal(1:end).*rhog./rhom.internal)+V0.*(1-alphag.internal(1:end)).*alphag.internal(1:end).*(rhog./rhom.internal-1)).*alphag.internal(1:end);

      flux.left.value=(V0*(1-alphag.left.setvalue)*(1-alphag.left.setvalue*rhog/rhom.left.setvalue)+V0*(1-alphag.left.setvalue)*alphag.left.setvalue.*(rhog/rhom.left.setvalue-1))*alphag.left.setvalue;
      
      flux.right.value=(V0*(1-alphag.right.setvalue)*(1-alphag.right.setvalue*rhog/rhom.right.setvalue)+V0*(1-alphag.right.setvalue)*alphag.right.setvalue.*(rhog/rhom.right.setvalue-1)).*alphag.right.setvalue;
      
      flux=setBC(flux,constField(0,N),xC,xF,0);
    
  else

      % rhom field
      rhom=assign(assign(alphag,constField(rhog,N),'*'),assign(assign(constField(1,N),alphag,'-'),constField(rhol,N),'*'),'+');
      rhom=setBC(rhom,constField(0,N),xC,xF,0);
  
      % Velocities actualization
      Vm.internal=V0.*(1-alphag.internal(1:end)).*alphag.internal(1:end).*(rhog./rhom.internal-1);  
      Vm=setBC(Vm,constField(0,N),xC,xF,0);

      % Calculation of drift velocity from alphag
      Vdrp.internal=V0.*(1-alphag.internal(1:end)).*(1-alphag.internal(1:end).*rhog./rhom.internal);  
      Vdrp=setBC(Vdrp,constField(0,N),xC,xF,0);  

      % Flux evaluation (remember that it is a cell flux, flux at boundaries are stored in BC's)
  
      flux=assign(assign(Vm,Vdrp,'+'),alphag,'*');
   
      flux=setBC(flux,constField(0,N),xC,xF,0);

  end

      cflux=zeros(N+1,1);

  if 0  
       

	  alphag=[alphagLeft; alphag; alphagRight];
	
	  eig1=-(sqrt((1-alphag).^2.*alphag.^2.*(1000000.*(15.*alphag.^6-20.*alphag.^5+6.*alphag.^4)+1000000000000.*(15.*alphag.^6-40.*alphag.^5+36.*alphag.^4-12.*alphag.^3+alphag.^2)+1000000000000000000.*(alphag.^6-4.*alphag.^5+6.*alphag.^4-4.*alphag.^3+alphag.^2)+alphag.^6+1000000000000000.*(-6.*alphag.^6+20.*alphag.^5-24.*alphag.^4+12.*alphag.^3-2.*alphag.^2)+1000.*(4.*alphag.^5-6.*alphag.^6)+1000000000.*(-20.*alphag.^6+40.*alphag.^5-24.*alphag.^4+4.*alphag.^3)+(alphag+1000.*(1-alphag)).*(1000000000.*(20.*alphag.^5-48.*alphag.^4+36.*alphag.^3-8.*alphag.^2)+1000.*(10.*alphag.^5-8.*alphag.^4)+1000000000000000.*(2.*alphag.^5-8.*alphag.^4+12.*alphag.^3-8.*alphag.^2+2.*alphag)-2.*alphag.^5+1000000000000.*(-10.*alphag.^5+32.*alphag.^4-36.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000.*(-20.*alphag.^5+32.*alphag.^4-12.*alphag.^3))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(6.*alphag.^4-12.*alphag.^3+6.*alphag.^2)+1000000000000.*(alphag.^4-4.*alphag.^3+6.*alphag.^2-4.*alphag+1)+alphag.^4+1000000000.*(-4.*alphag.^4+12.*alphag.^3-12.*alphag.^2+4.*alphag)+1000.*(4.*alphag.^3-4.*alphag.^4))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).^2.*V0.^2+0.282.*(1-alphag).*alphag.*((alphag+1000.*(1-alphag)).*(1000000.*(16.*alphag.^6-30.*alphag.^5+16.*alphag.^4-2.*alphag.^3)+1000000000000.*(16.*alphag.^6-58.*alphag.^5+80.*alphag.^4-52.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000000000000.*(-4.*alphag.^6+18.*alphag.^5-32.*alphag.^4+28.*alphag.^3-12.*alphag.^2+2.*alphag)+1000.*(4.*alphag.^5-4.*alphag.^6)+1000000000.*(-24.*alphag.^6+66.*alphag.^5-64.*alphag.^4+26.*alphag.^3-4.*alphag.^2))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(12.*alphag.^5-26.*alphag.^4+16.*alphag.^3-2.*alphag.^2)+1000000000000.*(4.*alphag.^5-18.*alphag.^4+32.*alphag.^3-28.*alphag.^2+12.*alphag-2)+1000.*(4.*alphag.^4-4.*alphag.^5)+1000000000.*(-12.*alphag.^5+40.*alphag.^4-48.*alphag.^3+24.*alphag.^2-4.*alphag))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.079524.*(alphag+1000.*(1-alphag)).^2.*(1000000.*(4.*alphag.^6-8.*alphag.^5+4.*alphag.^4)+1000000000000.*(4.*alphag.^6-20.*alphag.^5+41.*alphag.^4-44.*alphag.^3+26.*alphag.^2-8.*alphag+1)+1000000000.*(-8.*alphag.^6+28.*alphag.^5-36.*alphag.^4+20.*alphag.^3-4.*alphag.^2)))+(1-alphag).*alphag.*(1000000.*(3.*alphag.^3-4.*alphag.^2+alphag)+alphag.^3+1000000000.*(-alphag.^3+2.*alphag.^2-alphag)+1000.*(2.*alphag.^2-3.*alphag.^3)+(alphag+1000.*(1-alphag)).*(1000.*(6.*alphag.^2-6.*alphag)-3.*alphag.^2+1000000.*(-3.*alphag.^2+6.*alphag-3))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.282.*(alphag+1000.*(1-alphag)).*(1000000.*(2.*alphag.^3-5.*alphag.^2+4.*alphag-1)+1000.*(2.*alphag.^2-2.*alphag.^3)))./((alphag+1000.*(1-alphag)).*(1000000.*(2.*alphag.^2-4.*alphag+2)+2.*alphag.^2+1000.*(4.*alphag-4.*alphag.^2)));
	
	  eig2=(sqrt((1-alphag).^2.*alphag.^2.*(1000000.*(15.*alphag.^6-20.*alphag.^5+6.*alphag.^4)+1000000000000.*(15.*alphag.^6-40.*alphag.^5+36.*alphag.^4-12.*alphag.^3+alphag.^2)+1000000000000000000.*(alphag.^6-4.*alphag.^5+6.*alphag.^4-4.*alphag.^3+alphag.^2)+alphag.^6+1000000000000000.*(-6.*alphag.^6+20.*alphag.^5-24.*alphag.^4+12.*alphag.^3-2.*alphag.^2)+1000.*(4.*alphag.^5-6.*alphag.^6)+1000000000.*(-20.*alphag.^6+40.*alphag.^5-24.*alphag.^4+4.*alphag.^3)+(alphag+1000.*(1-alphag)).*(1000000000.*(20.*alphag.^5-48.*alphag.^4+36.*alphag.^3-8.*alphag.^2)+1000.*(10.*alphag.^5-8.*alphag.^4)+1000000000000000.*(2.*alphag.^5-8.*alphag.^4+12.*alphag.^3-8.*alphag.^2+2.*alphag)-2.*alphag.^5+1000000000000.*(-10.*alphag.^5+32.*alphag.^4-36.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000.*(-20.*alphag.^5+32.*alphag.^4-12.*alphag.^3))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(6.*alphag.^4-12.*alphag.^3+6.*alphag.^2)+1000000000000.*(alphag.^4-4.*alphag.^3+6.*alphag.^2-4.*alphag+1)+alphag.^4+1000000000.*(-4.*alphag.^4+12.*alphag.^3-12.*alphag.^2+4.*alphag)+1000.*(4.*alphag.^3-4.*alphag.^4))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).^2.*V0.^2+0.282.*(1-alphag).*alphag.*((alphag+1000.*(1-alphag)).*(1000000.*(16.*alphag.^6-30.*alphag.^5+16.*alphag.^4-2.*alphag.^3)+1000000000000.*(16.*alphag.^6-58.*alphag.^5+80.*alphag.^4-52.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000000000000.*(-4.*alphag.^6+18.*alphag.^5-32.*alphag.^4+28.*alphag.^3-12.*alphag.^2+2.*alphag)+1000.*(4.*alphag.^5-4.*alphag.^6)+1000000000.*(-24.*alphag.^6+66.*alphag.^5-64.*alphag.^4+26.*alphag.^3-4.*alphag.^2))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(12.*alphag.^5-26.*alphag.^4+16.*alphag.^3-2.*alphag.^2)+1000000000000.*(4.*alphag.^5-18.*alphag.^4+32.*alphag.^3-28.*alphag.^2+12.*alphag-2)+1000.*(4.*alphag.^4-4.*alphag.^5)+1000000000.*(-12.*alphag.^5+40.*alphag.^4-48.*alphag.^3+24.*alphag.^2-4.*alphag))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.079524.*(alphag+1000.*(1-alphag)).^2.*(1000000.*(4.*alphag.^6-8.*alphag.^5+4.*alphag.^4)+1000000000000.*(4.*alphag.^6-20.*alphag.^5+41.*alphag.^4-44.*alphag.^3+26.*alphag.^2-8.*alphag+1)+1000000000.*(-8.*alphag.^6+28.*alphag.^5-36.*alphag.^4+20.*alphag.^3-4.*alphag.^2)))+(1-alphag).*alphag.*(1000.*(3.*alphag.^3-2.*alphag.^2)+1000000000.*(alphag.^3-2.*alphag.^2+alphag)-alphag.^3+1000000.*(-3.*alphag.^3+4.*alphag.^2-alphag)+(alphag+1000.*(1-alphag)).*(1000000.*(3.*alphag.^2-6.*alphag+3)+3.*alphag.^2+1000.*(6.*alphag-6.*alphag.^2))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.282.*(alphag+1000.*(1-alphag)).*(1000.*(2.*alphag.^3-2.*alphag.^2)+1000000.*(-2.*alphag.^3+5.*alphag.^2-4.*alphag+1)))./((alphag+1000.*(1-alphag)).*(1000000.*(2.*alphag.^2-4.*alphag+2)+2.*alphag.^2+1000.*(4.*alphag-4.*alphag.^2)));
	  
	  % Local velocities at internal field
	  a=max(abs(eig1),abs(eig2));
	  
	  aLeft=a(1);
	  aRight=a(end);
	  % a is turned back to internal values only
	  a=a(2:end-1);
	  
	  % alphag is turned back to internal values only
	  alphag=alphag(2:end-1);
	
    
    a_j_minus_half=[aLeft;max(a(1:end-1),a(2:end))];
    a_j_plus_half=[max(a(1:end-1),a(2:end));aRight];
  
  else
  
    [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag,rhol,rhog,V0);
  
  end
  
  if 0
  
    % Rosunov scheme for non boundary cells
    alphag(2:end-1)=alphag(2:end-1)-lambda/2*(flux(3:end)-flux(1:end-2))+1/2*lambda*(a_j_plus_half(2:end-1).*(alphag(3:end)-alphag(2:end-1))-a_j_minus_half(2:end-1).*(alphag(2:end-1)-alphag(1:end-2)));
    
    % BC's
    alphag(1)=alphag(1)-lambda*(flux(2)-fluxleft)+lambda*(a_j_plus_half(1)*(alphag(2)-alphag(1))-a_j_minus_half(1)*(alphag(1)-alphagLeft));
    alphag(end)=alphag(end)-lambda*(fluxright-flux(end-1))+lambda*(a_j_plus_half(end)*(alphagRight-alphag(end))-a_j_minus_half(end)*(alphag(end)-alphag(end-1)));

  else
 
    [alphag,alphag]=Rusanov(alphag,alphag,flux,flux,cflux,cflux,a_j_minus_half,a_j_plus_half,0,0,dummyRho,dummyRho,dx,dt);	   
    %[alphag,alphag]=Rusanov(u1_0,u2_0,flux1,flux2,cflux1,cflux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,dx,dt)

  end

  % Calculation of mean velocity from alphag
  %Vm=V0.*(1-alphag.internal(1:end)).*alphag.internal(1:end).*(rhog./rhom.internal-1);  
  
  if 0
  
    figure(1); plot(xC,alphag);
    figure(2); plot(eig1);
    figure(3); plot(eig2);
    figure(4); plot(xC,Vm);
  
  end

end
