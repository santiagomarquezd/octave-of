function [eigens]=aspeed(U,alphag,rhol,rhog,N)
    % Gives the max of |eigenvalues| at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem (no Vr, no dp/dx)
    %
    % dU/dt+d(rhom*U*U)/dx=0
    % d alphag/dx+d(alphag*U)/dx=0
    %
    % [eigens]=a(U,alphag)
    %
    % eigens: max eigenvalues at interfaces
    % U: mixture velocity
    % rhol: liquid density
    % rhog: gas density
    % alphag: gas volume fraction
    
    % Allocation
    % The eigenvalues are defined at cell interfaces
    eigens=zeros(N+1,1);
    
    % Internal field caculation
    eigens(2:end-1)=max(max(abs(-(sqrt(4.*alphag.internal(2:end).*rhol-4.*alphag.internal(2:end).*rhog+1).*U.internal(2:end)-3.*U.internal(2:end))./2),abs((sqrt(4.*alphag.internal(2:end).*rhol-4.*alphag.internal(2:end).*rhog+1).*U.internal(2:end)+3.*U.internal(2:end))./2)),max(abs(-(sqrt(4.*alphag.internal(1:end-1).*rhol-4.*alphag.internal(1:end-1).*rhog+1).*U.internal(1:end-1)-3.*U.internal(1:end-1))./2),abs((sqrt(4.*alphag.internal(1:end-1).*rhol-4.*alphag.internal(1:end-1).*rhog+1).*U.internal(1:end-1)+3.*U.internal(1:end-1))./2)));
    
    % BC's
    eigens(1)=max(abs(-(sqrt(4*alphag.left.setvalue*rhol-4*alphag.left.setvalue*rhog+1)*U.left.setvalue-3*U.left.setvalue)/2),abs((sqrt(4*alphag.left.setvalue*rhol-4*alphag.left.setvalue*rhog+1)*U.left.setvalue+3*U.left.setvalue)/2));
    
    eigens(end)=max(abs(-(sqrt(4*alphag.right.setvalue*rhol-4*alphag.right.setvalue*rhog+1)*U.right.setvalue-3*U.right.setvalue)/2),abs((sqrt(4*alphag.right.setvalue*rhol-4*alphag.right.setvalue*rhog+1)*U.right.setvalue+3*U.right.setvalue)/2));
    
end

