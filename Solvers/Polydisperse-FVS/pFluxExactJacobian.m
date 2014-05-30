function [J]=pFluxExactJacobian(ui,uipo,V0,a)
    % Gives exact Jacobian for a system of nonlinear advection
    % equations in polydisperse sedimentation
    % The jacobian is evaluated in the ith,(i+1)th intercell
    % and stored in i position
    %
    % [J]=pFluxJacobian(ui,uipo,V0,a)
    %
    % J: calculated jacobian
    % ui: states vector of dispersed phases at cell i
    % uipo: states vector of dispersed phases at cell i+1
    % V0: relative velocity law constant
    % a: exponents vector

    % Number of dispersed phases
    N=size(ui,1);   
    
    % Memory allocation
    J=zeros(2,2);
 
    umean=(ui(1,1)+uipo(1,1))/2;
    J11=(1-umean)^2*V0(1,1)-2*(1-umean)*umean*V0(1,1);
    %J11=(1-umean)*V0(1,1)-umean*V0(1,1);
    umean=(ui(2,1)+uipo(2,1))/2;
    J22=(1-umean)^2*V0(2,1)-2*(1-umean)*umean*V0(2,1);
    %J22=(1-umean)*V0(2,1)-umean*V0(2,1);
    J=[J11 0; 0 J22];
end

    