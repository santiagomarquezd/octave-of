function [J]=arrayPFluxJacobian(u,V0,a,alphaDPL)
    % Gives numerical Jacobian  for a system of nonlinear advection
    % equations in polydisperse sedimentation in array form
    % (forward differencing)
    %
    % [J]=pFluxJacobian(u,V0,a)
    %
    % J: array of calculated jacobians
    % u: states vector array of dispersed phases
    % V0: relative velocity law constant
    % a: exponents vector
    % alphaDPL: concentration of the dense-packed layer

    % Number of dispersed phases and points for jacobian calculation
    [N,M]=size(u);   
    
    % Memory allocation
    J=zeros(N,N,M);

    % u variable lives in [0, 1]
    du=1/1000; 

    % Jacobian calculation loop
    % Derivative of fluxes 1:N respect to variable i
    if (1)
        % Forward differencing
        for i=1:N
            Fu=arrayPFlux(u,V0,a,alphaDPL);
            uPlusdu=u;
            uPlusdu(i,:)=u(i,:)+du;
            FuPlusdu=arrayPFlux(uPlusdu,V0,a,alphaDPL);
            J(:,i,:)=(FuPlusdu-Fu)/du;
        end
    else
        % Central differencing
        for i=1:N
            uMinusdu=u;
            uMinusdu(i,:)=u(i,:)-du/2;
            FuMinusdu=arrayPFlux(uMinusdu,V0,a,alphaDPL);
            uPlusdu=u;
            uPlusdu(i,:)=u(i,:)+du/2;
            FuPlusdu=arrayPFlux(uPlusdu,V0,a,alphaDPL);
            J(:,i,:)=(FuPlusdu-FuMinusdu)/du;
        end
    end
end
