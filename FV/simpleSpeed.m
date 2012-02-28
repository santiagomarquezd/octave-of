function [eigens]=simpleSpeed(V0, alphag,rhol,rhog, rhom,N)
    % Gives the max eigenvalues at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem (algebraic U)
    %
    % d alphag/dx+d(alphag*Vg)/dx=0
    %
    % [eigens]=simpleSpeed(V0,alphag,rhol,rhog,N)
    %
    % eigens: max eigenvalues at interfaces
    % V0: multiplicative constant in Vpq equation
    % alphag: gas volume fraction
    % rhol: liquid density
    % rhog: gas density
    % rhom: mixture density
    % N: number of cells
    
    % Allocation
    % The eigenvalues are defined at cell interfaces
    eigens=zeros(N+1,1);
    
    % Internal field caculation
    eigens(2:end-1)=max(abs(-1+2.*alphag.internal(2:end)-(1-alphag.internal(2:end).*rhog./rhom.internal(2:end))+(1-alphag.internal(2:end)).*(-rhog.*rhom.internal(2:end)+alphag.internal(2:end).*rhog.*(rhog-rhol))./rhom.internal(2:end).^2),abs(abs(-1+2.*alphag.internal(2:end)-(1-alphag.internal(1:end-1).*rhog./rhom.internal(1:end-1))+(1-alphag.internal(1:end-1)).*(-rhog.*rhom.internal(1:end-1)+alphag.internal(1:end-1).*rhog.*(rhog-rhol))./rhom.internal(1:end-1).^2);
    
    % BC's
    eigens(1)=abs(-1+2.*alphag.left.setvalue-(1-alphag.left.setvalue.*rhog./rhom.left.setvalue)+(1-alphag.left.setvalue).*(-rhog.*rhom.left.setvalue+alphag.left.setvalue.*rhog.*(rhog-rhol))./alphag.left.setvalue.^2));
    
    eigens(end)=abs(-1+2.*alphag.right.setvalue-(1-alphag.right.setvalue.*rhog./rhom.right.setvalue)+(1-alphag.right.setvalue).*(-rhog.*rhom.right.setvalue+alphag.right.setvalue.*rhog.*(rhog-rhol))./alphag.right.setvalue.^2));
    
end

