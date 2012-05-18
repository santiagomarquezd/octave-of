function [u1,u2]=FVS(u1_0,u2_0,Aplus,Aminus,dx,dt,N,imp)
  % Applies the Flux Vector Splitting Method to U=[u1;u2]
  % using (8.47) and (8.48) of Riemann Solvers and Numerical Methods
  % for Fluid Dynamics (2nd. edition) from E.F.Toro
  %
  %
  % [u1,u2]=FVS(u1_0,u2_0,Aplus,Aminus,dx,dt,imp)
  %
  % u1: u1 scalar unknown time advanced
  % u2: u2 scalar unknown time advanced
  % u1_0: u1 scalar unknown at previous time-step
  % u2_0: u2 scalar unknown at previous time-step
  % Aplus: array of positive eigenvalues matrices
  % Aminus: array of negative eigenvalues matrices
  % dx: mesh step
  % dt: time-step
  % N: number of cells
  % imp:indicates wheter the walls are impermeable
  %     1, impermeable

  % Memory allocation
  u1=u1_0;
  u2=u2_0;

  % Data reshaping
  % States as column vectors (internal values)
  U0=[u1_0.internal';u2_0.internal'];
  
  % Memory allocation
  % States
  U=U0;
  % Fluxes
  F=zeros(2,N-1);

  % Splitted fluxes calculation
  % Internal faces
  for i=1:(N-1)
    F(:,i)=Aplus(:,:,i)*U0(:,i)+Aminus(:,:,i+1)*U0(:,i+1);
  end
    
  % Boundary conditions
  Fleft=Aplus(:,:,i)*[u1_0.left.setvalue;u2_0.left.setvalue]+Aminus(:,:,1)*U0(:,1);
  Fright=Aminus(:,:,N)*[u1_0.right.setvalue;u2_0.right.setvalue]+Aplus(:,:,N)*U0(:,N);

  % Sets impermeable walls
  if (imp)
    Fleft=zeros(2,1);
    Fright=zeros(2,1);
  end

  % Final flux	
  F=[Fleft F Fright];

  

  % Time advancement
  for i=1:N
    U(:,i)=U(:,i)-dt/dx*(F(:,i+1)-F(:,i));
  end

  %keyboard; pause;

  % Data reshaping to external format
  u1.internal=U(1,:)';
  u2.internal=U(2,:)';

end