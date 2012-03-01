function [tddtCouplingCoeff] = fvc_DdtPhiCoeff_3arg(U, phi, phiCorr)
    % Gives the coefficient for fvcDdtPhi correction calculation
    %
    % [tddtCouplingCoeff] = fvc_DdtPhiCoeff_3arg(U, phi, phiCorr)
    %
    % U: velocity field
    % phi: volumetric flux at faces
    % phiCorr: volumetric flux corrected at faces
    
% ddtScheme<Type>::fvcDdtPhiCoeff
%  (
%      const GeometricField<Type, fvPatchField, volMesh>& U,
%      const fluxFieldType& phi,
%      const fluxFieldType& phiCorr
%  )
%  {
%      tmp<surfaceScalarField> tddtCouplingCoeff = scalar(1)
%        - min
%          (
%              mag(phiCorr)
%             /(mag(phi) + dimensionedScalar("small", phi.dimensions(), SMALL)),
%              scalar(1)
%          );
%  
%      surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff();
%  
%      forAll(U.boundaryField(), patchi)
%      {
%          if (U.boundaryField()[patchi].fixesValue())
%          {
%              ddtCouplingCoeff.boundaryField()[patchi] = 0.0;
%          }
%      }
%  
%      if (debug > 1)
%      {
%          Info<< "ddtCouplingCoeff mean max min = "
%              << gAverage(ddtCouplingCoeff.internalField())
%              << " " << gMax(ddtCouplingCoeff.internalField())
%              << " " << gMin(ddtCouplingCoeff.internalField())
%              << endl;
%      }
%  
%      return tddtCouplingCoeff;
%  }

tddtCouplingCoeff = 1 -min(...
                           abs(phiCorr)/(abs(phi)+ 1E-37)...
                           ,1);
                
 
    if (U.left.type=='V')
      tddtCouplingCoeff(1)=0;
    elseif (U.right.type=='V')
      tddtCouplingCoeff(end)=0;
    end