function [fluxFlat]=AEqnNoVmFluxFlat(uFlat,V0,rhol,rhog,a)
    % Gives the flux for alphaEqn written in A in mixture model
    % with no Vm term u flatten version (u has plain data)
    % IMPORTANT!: the flux is calculated using alphag
    %
    % [fluxFlat]=AEqnNoVmFluxFlat(u,Vpq,rhom,cp)
    %
    % fluxFlat: flatten obtained flux
    % u: independent unknown in flat form
    % V0: constant for relative velocity calculation
    % rhol: continue phase density 
    % rhog: dispersed phase density
    % a: exponent for relative velocity calculation 

    % Deflatten
    u.internal=uFlat;
    u.left.setvalue=uFlat(1);
    u.right.setvalue=uFlat(end);

    % Takes the data size
    N=size(u.internal,1);
    
    rhom=assign(assign(constField(rhog,N),u,'*'),assign(assign(constField(1,N),u,'-'),constField(rhol,N),'*'),'+');
    cp=assign(assign(constField(rhog,N),u,'*'),rhom,'/');

    % Flux=V0*(1-alphag)^a*(1-cp)*rhom*A=V0*(1-alphag)^a*(1-cp)*rhog*alphag
    flux=assign(assign(assign(assign(constField(V0,N),assign(assign(constField(1,N),u,'-'),constField(a,N),'^'),'*'),assign(constField(1,N),cp,'-'),'*'),constField(rhog,N),'*'),u,'*');
    % The flux results to be rhom, rhol and rhog independent: ((1-alphag)^a*alphag-(1-alphag)^a*alphag^2)*V0
    %flux=assign(assign(assign(assign(assign(constField(1,N),u,'-'),constField(a,N),'^'),u,'*'),assign(assign(assign(constField(1,N),u,'-'),constField(a,N),'^'),assign(u,constField(2,N),'^'),'*'),'-'),constField(V0,N),'*');
 
    % Flatten
    fluxFlat=flux.internal;
end