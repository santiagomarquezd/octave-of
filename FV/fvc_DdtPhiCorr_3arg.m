function [fvcDdtPhiCorr] = fvc_DdtPhiCorr_3arg(rA, U, phiAbs,w,xC,xF,Sf,dt,N)
    %    % Explicit correction of mass flux taken into account changes in time
    %
    % [fvcDdtPhiCorr] = fvc_DdtPhiCorr_3arg(rA, U, phiAbs,w,xC,xF,Sf,dt,N)
    %
    % rUA: 1/Aoperator
    % U: velocity at previous timestep
    % phiabs: volumetric flux at previous timestep
    % w: interpolation weights
    % xC: cell centres positions
    % xF: face centres positions
    % Sf: face areas
    % dt: timestep  


% EulerDdtScheme<Type>::fvcDdtPhiCorr
% (
%     const volScalarField& rA,
%     const GeometricField<Type, fvPatchField, volMesh>& U,
%     const fluxFieldType& phiAbs
% )
% {
%     dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
rDeltaT = 1/dt;

%     IOobject ddtIOobject
%     (
%         "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phiAbs.name() + ')',
%         mesh().time().timeName(),
%         mesh()
%     );
% 
%     tmp<fluxFieldType> phiCorr =
%         phiAbs.oldTime() - (fvc::interpolate(U.oldTime()) & mesh().Sf());

phiCorr = phiAbs - (fvc_interpolate(U, w, xC, xF).*Sf);

%     return tmp<fluxFieldType>
%     (
%         new fluxFieldType
%         (
%             ddtIOobject,
%             this->fvcDdtPhiCoeff(U.oldTime(), phiAbs.oldTime(), phiCorr())
%           * fvc::interpolate(rDeltaT*rA)*phiCorr
%         )
%     );
% }

fvcDdtPhiCorr = fvc_DdtPhiCoeff_3arg(U, phiAbs,phiCorr)* fvc_interpolate(assign(constField(rDeltaT,N),rA,'*'),w).* phiCorr;
