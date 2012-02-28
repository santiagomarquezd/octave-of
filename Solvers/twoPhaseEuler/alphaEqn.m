% alphaEqn
%  {
%      word scheme("div(phi,alpha)");
%      word schemer("div(phir,alpha)");
%  
%      surfaceScalarField phic("phic", phi);
%      surfaceScalarField phir("phir", phia - phib);
%
phic = phi;
phir = phia0 - phib0;

%      if (g0.value() > 0.0)
%      {
%          surfaceScalarField alphaf(fvc::interpolate(alpha));
%          surfaceScalarField phipp(ppMagf*fvc::snGrad(alpha)*mesh.magSf());
%          phir += phipp;
%          phic += fvc::interpolate(alpha)*phipp;
%      }
%  
%      for (int acorr=0; acorr<nAlphaCorr; acorr++)
%      {
%          fvScalarMatrix alphaEqn
%          (
%               fvm::ddt(alpha)
%             + fvm::div(phic, alpha, scheme)
%             + fvm::div(-fvc::flux(-phir, beta, schemer), alpha, schemer)
%          );
%  
%          if (g0.value() > 0.0)
%          {
%              ppMagf = rUaAf*fvc::interpolate
%              (
%                  (1.0/(rhoa*(alpha + scalar(0.0001))))
%                 *g0*min(exp(preAlphaExp*(alpha - alphaMax)), expMax)
%              );
%  
%              alphaEqn -= fvm::laplacian
%              (
%                  (fvc::interpolate(alpha) + scalar(0.0001))*ppMagf,
%                  alpha,
%                  "laplacian(alphaPpMag,alpha)"
%              );
%          }
%  
%          alphaEqn.relax();
%          alphaEqn.solve();
%  
%          #include "packingLimiter.H"
%  
%          beta = scalar(1) - alpha;
%  
%          Info<< "Dispersed phase volume fraction = "
%              << alpha.weightedAverage(mesh.V()).value()
%              << "  Min(alpha) = " << min(alpha).value()
%              << "  Max(alpha) = " << max(alpha).value()
%              << endl;
%      }
%  }
%  
%  rho = alpha*rhoa + beta*rhob;

% Assembling
disp('Assembling alphaGEqn')
[M, RHS] = fvm_ddt(constField(1,N),constField(1,N),alpha0,V,dt,1);
[A, ARHS]= fvm_div_flux_cell(phic,alpha0,xC,xF,w,1);
[B, BRHS]= fvm_div_flux_cell(phir.*fvc_interpolate(beta0, w, xC, xF),alpha0,xC,xF,w,1);  

M=M+A+B;
RHS=RHS+ARHS+BRHS;
 
% Solve
disp('Solving for alphag')
alpha.internal=M\RHS;

beta = assign(constField(1,N),alpha,'-');
rho = assign(assign(alpha,constField(transportProperties.phasea.rho,N),'*'),assign(beta,constField(transportProperties.phaseb.rho,N),'*'),'+');
