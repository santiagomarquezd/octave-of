function [fluxMax,fluxMin]=fluxMaxMin(u1,u2,fluxFunc,rhol,rhog,V0,a,N)
    % Gives the flux maximum and minimum for a given
    % interval of independent unknown
    %
    % [fluxMax,fluxMin]=fluxMaxMin(u1,u2,flux,rhog,rhol,V0)
    %
    % fluxMax: maximum of the flux in the interval
    % fluxMin: minimum of the flux in the interval
    % u1: first interval extremum
    % u1: second interval extremum
    % fluxFunc: flux function
    % rhol: liquid density
    % rhog: gas density
    % V0: constant for flux equation
    % a: exponent for flux equation
    % N: number of subinterval for extrama searching 
    %    bigger N gives more precise extrema
    
    % Generates a valid interval extrema
    if (u1<u2)
      uL=u1;
      uR=u2;
    elseif (u2<=u1)
      uL=u2;
      uR=u1;
    end

    % Generates interval
    if (u1==u2)
      % In case of equal values the interval are those values
      I=[u1 u2];
    else
      du=(uR-uL)/N;
      I=uL:du:uR;
    end  

    % Field generation
    u.internal=I';
    % Dummy boundary conditions, only internal field is needed
    u.left.type='V';
    u.left.setvalue=0;
    u.right.type='V';
    u.right.setvalue=0;

    %keyboard; pause;

    % Flux calculation
    F=fluxFunc(u,V0,rhol,rhog,a);    

    % Finds maximum and minimum
    fluxMax=max(F.internal);
    fluxMin=min(F.internal);    

end