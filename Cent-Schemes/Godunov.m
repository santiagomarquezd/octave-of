function [u1]=Godunov(u1_0,flux1,cflux1,rhom0,rhom,V0,rhol,rhog,a,dx,dt,Nint)
  % Applies the Godunov scheme given actual state and cell fluxes
  %
  % [u1,u2]=Rusanov(u1_0,u2_0,flux1,flux2,cflux1,cflux2,a_j_minus_half,a_j_plus_half,S1,S2,rhom0,rhom,dx,dt)
  %

  % Data size
  N=size(u1_0.internal,1);

  % Allocation
  u1=u1_0;
  F=zeros(N+1,1);

  % Calculates faces flux by Godunov's method
  for i=1:N-1
    [fluxMax,fluxMin]=fluxMaxMin(u1_0.internal(i,1),u1_0.internal(i+1,1),flux1,rhol,rhog,V0,a,Nint);
    if (u1_0.internal(i,1)<=u1_0.internal(i+1,1))
      F(i+1)=fluxMin;
    elseif (u1_0.internal(i,1)>=u1_0.internal(i+1,1))
      F(i+1)=fluxMax;
    end
  end
   
  % Impermeable walls
  F(1)=0;
  F(end)=0;

  % Explicit temporal scheme
  u1.internal(1:end)=u1_0.internal(1:end).*rhom0.internal./rhom.internal-dt./dx*(F(2:end)-F(1:end-1))./rhom.internal;

end
  
