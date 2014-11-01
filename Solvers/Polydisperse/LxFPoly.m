function [u1,u2,u3]=LxFPoly(u1_0,u2_0,u3_0,V0,a,alphaDPL,dx,dt)
    % Applies the Lax-Friederichs scheme for polydisperse cases
    %
    % [u1,u2]=LxFPoly(u1_0,u2_0,flux1,flux2,dx,dt)
    %
    % u1,u2,u3: new states (only internal field)
    % u1_0,u2_0,u3_0: previous states (only internal field)
    % V0: relative velocity law constant
    % a: exponents vector
    % alphaDPL: concentration of the dense-packed layer
    % dx: spatial step
    % dt: time-step

    % Get cells number
    N=size(u1_0,1);

    % Data allocation
    u1=u1_0;
    u2=u2_0;
    u3=u3_0;
    
    % Temporal data
    u1tmp=u1_0;
    u2tmp=u2_0;
    u3tmp=u3_0;
    usum=u1_0+u2_0+u3_0;

    % Data reshaping
    u=[ (u1tmp)' 
        (u2tmp)' 
        (u3tmp)'];
    
    F=arrayPFlux(u,V0,a,alphaDPL);

    F1=F(1,:)';
    F2=F(2,:)';
    F3=F(3,:)';
    
    % u1 temporal advancement
    % Non boundary cells
    u1(2:N-1)=1/2*(u1tmp(3:N)+u1tmp(1:N-2))-dt/dx/2*(F1(3:N)-F1(1:N-2));

    % Boundary cells           
    % First cell
    u1(1)=1/2*(u1tmp(1)+u1tmp(2))-dt/dx/2*(F1(1)+F1(2)); % Impermeable wall
    % Last cell
    u1(N)=1/2*(u1tmp(N-1)+u1tmp(N))+dt/dx/2*(F1(N-1)+F1(N)); % Impermeable wall


    % u2 temporal advancement
    % Non boundary cells
    u2(2:N-1)=1/2*(u2tmp(3:N)+u2tmp(1:N-2))-dt/dx/2*(F2(3:N)-F2(1:N-2));

    % Boundary cells           
    % First cell
    u2(1)=1/2*(u2tmp(1)+u2tmp(2))-dt/dx/2*(F2(1)+F2(2)); % Impermeable wall
    % Last cell
    u2(N)=1/2*(u2tmp(N-1)+u2tmp(N))+dt/dx/2*(F2(N-1)+F2(N)); % Impermeable wall

    % u3 temporal advancement
    % Non boundary cells
    u3(2:N-1)=1/2*(u3tmp(3:N)+u3tmp(1:N-2))-dt/dx/2*(F3(3:N)-F3(1:N-2));

    % Boundary cells           
    % First cell
    u3(1)=1/2*(u3tmp(1)+u3tmp(2))-dt/dx/2*(F3(1)+F3(2)); % Impermeable wall
    % Last cell
    u3(N)=1/2*(u3tmp(N-1)+u3tmp(N))+dt/dx/2*(F3(N-1)+F3(N)); % Impermeable wall

end