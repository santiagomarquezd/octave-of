% volVectorField Ur = Ua - Ub;
% volScalarField magUr = mag(Ur);

function [K] = dragmodel_K(Ur,InterfacialProperties,transportProperties)

 if (max(InterfacialProperties.dragModela == 'SchillerNaumann'))
    da = transportProperties.phasea.d;
    nub = transportProperties.phaseb.nu;
    rhob = transportProperties.phaseb.rho;
    
    magUr = Ur.internal;
    
    Re = max(magUr*da/nub, 1.0e-3);
    
  K = zeros(size(magUr,1),1);
  
  for i=1:size(magUr,1)
    if Re(i) <= 1000
        Cds = 24*(1 + 0.15*Re(i)^0.687)/Re(i);
    elseif Re(i) > 1000
        Cds = 0.44;
    end
    
    K(i) = 0.75*Cds*rhob*magUr(i)/da;
  end
 end
 
end
