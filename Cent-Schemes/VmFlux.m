function [flux]=VmFlux(Vm,alphag,Vr,rhol,rhog,model)
    % Gives the flux for Vm for several mixtureModels
    %
    % [flux]=VmFlux(alphag,Vm,Vr,rhol,rhog,model)
    %
    % flux: obtained flux
    % Vm: mixture velocity
    % alphag: alphag field
    % rhol: liquid density
    % rhog: gas density
    % model: 1. No Vr; 2: with Vr
    
    if (model==1)
	  rhom=alphag.*rhog+(1-alphag).*rhol;    
	  flux=rhom.*Vm.*Vm;      
    elseif (model==2)
	  disp('Not implemented yet')
    end
end