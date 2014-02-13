% calcVdrp
if (VpqModel == 1)

  % Implements relative velocity model from UADE    
   Vpq.internal = V0.*((alphaMax-min(alphag0.internal,alphaMax))/alphaMax).^aexp;
%    Vpq.internal=ones(N,1)*V0;
%    for i=1:N
%      if alphag.internal(i)>0.995
%        Vpq.internal(i) = 0;
%      end
%    end

elseif (VpqModel == 2)

  % Implements Eqn. 23.4-11 from Fluent's User's Guide
  % Precomputes needed values

  % Acceleration 23.4-15 from Fluent's User's Guide 
  a=g-U.internal.*fvc_grad(U, w, xC, xF, Sf, V)-fvc_ddt(U, U0, dt); 

  % Particle relaxation time 23.4-13 from Fluent's User's Guide 
  taup=rhog*dp^2/(18*mul);

  % Estimates Reynolds number using previous step relative velocity
  % (initialized as zero) Implements Manninen (39)
  % BC overriden.
  Re=dp*rhol*abs(Vpq.internal)/mul;

  % Schiller-Naumann's drag force 23.4-14 from Fluent's User's Guide 
  fdrag=Re*0;
  for i=1:size(Re,1)
    if (Re(i)<=1000)
      fdrag(i,1)=1+0.15*Re(i)^0.687;
    else
      fdrag(i,1)=0.0183*Re(i);
    end
  end

  % Relative velocity calculation Eqn. 23.4-12 from Fluent's User's Guide
  if 0
    Vpq.internal=taup./(fdrag+1E-6).*(rhog-rhom.internal)./rhog.*a;
  else
    Vpq.internal=-0.000037*(rhog-rhom.internal)./rhog*100;
  end

elseif (VpqModel == 3)

  Vpq.internal = V0*ones(N,1);

end

% To adjust BC's to dictionary
Vpq=setBC(Vpq,rhom,xC,xF,g);

% Mass fraction 23.4-10 from Fluent's User's Guide 
cp=alphag.internal.*rhog./rhom0.internal;

% Finally Drift Velocity calculation Eqn. 23.4-11 from Fluent's User's Guide
% using Eqn. 28 from Manninen
Vdrp=(1-cp).*Vpq.internal; %------->>>>>>>>>> VERRRRRR  Vdrp==(1-cp).*Vpq;



