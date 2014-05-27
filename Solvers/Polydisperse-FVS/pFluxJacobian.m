function [J]=pFluxJacobian(u,V0,a)
    % Gives numerical Jacobian  for a system of nonlinear advection
    % equations in polydisperse sedimentation
    % (forward differencing)
    %
    % [J]=pFluxJacobian(u,V0,a)
    %
    % J: calculated jacobian
    % u: states vector of dispersed phases
    % V0: relative velocity law constant
    % a: exponents vector

    % Number of dispersed phases
    N=size(u,1);   
    
    % Memory allocation
    J=zeros(N,N);

    % u variable lives in [0, 1]
    du=1/1000; 

    % Jacobian calculation loop
    for i=1:N
        % Derivative of fluxes 1:N respect to 
        % variable i
        Fu=pFlux(u,V0,a);
        uPlusdu=u;
        uPlusdu(i,1)=u(i,1)+du;
        FuPlusdu=pFlux(uPlusdu,V0,a);
        J(:,i)=(FuPlusdu-Fu)/du;
    end
end
