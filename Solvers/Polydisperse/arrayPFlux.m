function [F]=arrayPFlux(u,V0,a,alphaDPL)
    % Calculates the flux vector for a system of nonlinear advection
    % equations in polydisperse sedimentation (three classes)
    %
    % [F]=arrayPFlux(u,V0,a)
    %
    % F: calculated flux
    % u: array of states for dispersed phases
    % V0: relative velocity law constant
    % a: exponents vector
    % alphaDPL: concentration of the dense-packed layer

    % Array allocation
    F=u*0;

    % Calculation of beta
    usum=sum(u);

    u1=u(1,:);
    u2=u(2,:);
    u3=u(3,:);    

    F(1,:)=(V0(1,1).*(1-u1).*(alphaDPL-usum).^a(1,1)-1*(V0(2,1).*(alphaDPL-usum).^a(2,1).*u2+V0(3,1).*(alphaDPL-usum).^a(3,1).*u3)).*u1;
    F(2,:)=(V0(2,1).*(1-u2).*(alphaDPL-usum).^a(2,1)-1*(V0(1,1).*(alphaDPL-usum).^a(1,1).*u1+V0(3,1).*(alphaDPL-usum).^a(3,1).*u3)).*u2;
    F(3,:)=(V0(3,1).*(1-u3).*(alphaDPL-usum).^a(3,1)-1*(V0(1,1).*(alphaDPL-usum).^a(1,1).*u1+V0(2,1).*(alphaDPL-usum).^a(2,1).*u2)).*u3;
    
end

