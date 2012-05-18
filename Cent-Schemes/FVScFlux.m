function [u2]=FVScFlux(u1_0,u2_0,AAbs,Flux,cFlux,rhom,rhom0,dx,dt,N,xC,xF,w)
  % Applies the Flux Vector Splitting Method to u2 from U=[u1;u2]
  % using (5.76) and (5.105) of Riemann Solvers and Numerical Methods
  % for Fluid Dynamics (2nd. edition) from E.F.Toro adding
  % an external conservative flux. Only the second variable is updated
  % but the calculation depends in both u1 and u2
  %
  % [u2]=FVScFlux(u1_0,u2_0,AAbs,Flux,cFlux,rhom,rhom0,dx,dt,N,weights)
  %
  % u2: u2 scalar unknown time advanced
  % u1_0: u1 scalar unknown at previous time-step
  % AAbs: array of absolute value matrices
  % Flux: cell centered flux 
  % cFlux: conservative flux at faces 
  % rhom: density at advanced time (stimated)
  % rhom0: density at previous time
  % dx: mesh step
  % dt: time-step
  % N: number of cells
  % xC: cell centroids
  % xF: face centroids
  % w: weights for linear interpolation

  % Memory allocation
  u2=u1_0;

  % Data reshaping
  % States as column vectors (internal values)
  U0=[u1_0.internal';u2_0.internal'];
  % Densities 
  dummyRho=ones(N,1);
  Rho=[dummyRho';rhom0.internal'];

  %keyboard; pause;

  % Memory allocation
  % States
  U=U0;

  % Centered flux in faces from cell flux
  Fc=cFlux+fvc_interpolate(Flux, w, xC, xF);

  % Stabilization (only for internal faces)
  % values are multiplied by rhom0 for units consistency
  % when only the second unknown is multiplied
  Fs=zeros(N-1,1);
  for i=1:(N-1)
      Stab=-1/2*AAbs(:,:,i)*(U0(:,i+1).*Rho(:,i+1)-U0(:,i).*Rho(:,i));
      Fs(i)=Stab(2);
  end
 

  % Final flux
  F=[Fc(1);Fc(2:end-1)+Fs;Fc(end)];
  
  %keyboard; pause;

  % Time advancement
  u2.internal(1:end)=u2_0.internal(1:end).*rhom0.internal(1:end)./rhom.internal(1:end)-...
		     dt./dx*(F(2:end)-F(1:end-1))./rhom.internal(1:end);

end