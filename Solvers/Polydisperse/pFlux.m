function [F]=pFlux(u,V0,a)
    % Calculates the flux vector for a system of nonlinear advection
    % equations in polydisperse sedimentation
    %
    % [F]=FVSpFlux(u,a)
    %
    % F: calculated flux
    % u: states vector of dispersed phases
    % V0: relative velocity law constant
    % a: exponents vector

    % Number of dispersed phases
    N=size(u,1);   
    
    %keyboard; pause;

    % Calculation of beta
    beta=sum(u);

    % Fluxes calculation loop
    for i=1:N
        summ=0;
        for j=1:N
            if (j!=i)
                summ += u(j,1)*(V0(j,1)*(1-beta)^a(j,1));        
            end
        end
        %F(i,1) = u(i,1)*(V0(i,1)*(1-beta)^a(i,1)*(1-u(i,1))-summ);
        % Constant velocity flux
        %F(i,1) = u(i,1)*(V0(i,1));
        % Non constant velocity flux
        F(i,1) = u(i,1)*(V0(i,1))*(1-u(i,1))^a*(1-u(i,1));        
    end
end

