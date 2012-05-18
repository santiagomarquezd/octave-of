function [flux]=quadraticTrafficFlux(u,V0,rhol,rhog,a)
    % Gives the quadratic flux for Traffic Flow problem
    %
    % [flux]=quadraticTrafficFlux(u,V0)
    %
    % flux: obtained flux
    % u: independent unknown
    % V0: constant for flux calculation
    % rhol,rhog,a: for prototype compatibility
  
    % Takes the data size
    N=size(u.internal,1);

    %keyboard; pause;

    % Flux calculation    
    flux=assign(constField(V0,N),assign(assign(constField(1,N),u,'-'),u,'*'),'*');    
end