% UEqn.m

% fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
% fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);
% 
% {
%     {
%         volTensorField gradUaT(T(fvc::grad(Ua)));
% 
%         if (kineticTheory.on())
%         {
%             kineticTheory.solve(gradUaT);
%             nuEffa = kineticTheory.mua()/rhoa;
%         }
%         else // If not using kinetic theory is using Ct model
%         {
%             nuEffa = sqr(Ct)*nutb + nua;
%         }
nuEffa = transportProperties.phasea.nu*ones(N+1,1);
nuEffb = transportProperties.phasea.nu*ones(N+1,1);

%         volTensorField Rca
%         (
%             "Rca",
%             ((2.0/3.0)*I)*(sqr(Ct)*k + nuEffa*tr(gradUaT)) - nuEffa*gradUaT
%         );
% 
%         if (kineticTheory.on())
%         {
%             Rca -= ((kineticTheory.lambda()/rhoa)*tr(gradUaT))*tensor(I);
%         }
% 
%         surfaceScalarField phiRa
%         (
%             -fvc::interpolate(nuEffa)*mesh.magSf()*fvc::snGrad(alpha)
%             /fvc::interpolate(alpha + scalar(0.001))
%         );
%phiRa = -fvc_interpolate(nuEffa, w, xC, xF).*abs(Sf).*fvc_snGrad(alpha, w, xC, xF, Sf, V)/fvc_interpolate(alpha + 0.001, w, xC, xF);
phiRa = -nuEffa.*abs(Sf).*fvc_SnGrad(alpha, xC, xF)./fvc_interpolate(assign(alpha,constField(0.001,N),'+'), w, xC, xF);
% 
disp('Assembling UEqna')
%         UaEqn =
%         (
%             (scalar(1) + Cvm*rhob*beta/rhoa)*
%             (
%                 fvm::ddt(Ua)
%               + fvm::div(phia, Ua, "div(phia,Ua)")
%               - fvm::Sp(fvc::div(phia), Ua)
%             )
% 
%           - fvm::laplacian(nuEffa, Ua)
%           + fvc::div(Rca)
% 
%           + fvm::div(phiRa, Ua, "div(phia,Ua)")
%           - fvm::Sp(fvc::div(phiRa), Ua)
%           + (fvc::grad(alpha)/(fvc::average(alpha) + scalar(0.001)) & Rca)
%          ==
%         //  g                          // Buoyancy term transfered to p-equation
%           - fvm::Sp(beta/rhoa*K, Ua)
%         //+ beta/rhoa*K*Ub             // Explicit drag transfered to p-equation
%           - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)
%         );
% 
%         UaEqn.relax();
%     }
% 
[M, RHS] = fvm_ddt(constField(1,N),constField(1,N),Ua0,V,dt,1); % fvm::ddt(Ua)
[A, ARHS]= fvm_div_flux_cell(phia0,Ua0,xC,xF,w,1); % fvm::div(phia, Ua, "div(phia,Ua)")
ASp = fvm_Sp(fvc_div_face(phia0, V),V); % fvm::Sp(fvc::div(phia), Ua)

M = diag(ones(N,1)+(transportProperties.Cvm)*(transportProperties.phaseb.rho)*(beta0.internal)/(transportProperties.phasea.rho))*M;
A = diag(ones(N,1)+(transportProperties.Cvm)*(transportProperties.phaseb.rho)*(beta0.internal)/(transportProperties.phasea.rho))*A;
ASp = diag(ones(N,1)+(transportProperties.Cvm)*(transportProperties.phaseb.rho)*(beta0.internal)/(transportProperties.phasea.rho))*ASp;
RHS = (ones(N,1)+(transportProperties.Cvm)*(transportProperties.phaseb.rho)*(beta0.internal)/(transportProperties.phasea.rho)).*RHS;
ARHS = (ones(N,1)+(transportProperties.Cvm)*(transportProperties.phaseb.rho)*(beta0.internal)/(transportProperties.phasea.rho)).*ARHS;

[B, BRHS] = fvm_laplacian(Ua0,nuEffa,xC,xF,Sf); % fvm::laplacian(nuEffa, Ua)

[C, CRHS]= fvm_div_flux_cell(phiRa,Ua0,xC,xF,w,1); % fvm::div(phiRa, Ua, "div(phia,Ua)")
CSp = fvm_Sp(fvc_div_face(phiRa, V),V);           % fvm::Sp(fvc::div(phiRa), Ua)

%-------------- = -----------------------

% - fvm::Sp(beta/rhoa*dragCoef, Ua)
RHS1 = fvm_Sp(beta0.internal/(transportProperties.phasea.rho).*K, V);  % fvm::Sp(beta/rhoa*K, Ua)
RHS2 = beta.internal/(transportProperties.phaseb.rho).*(liftCoef + (transportProperties.Cvm)*(transportProperties.phaseb.rho)*DDtUb);

% Pressure terms are calculated apart in order to correctly calculate the H operator
%FRHS=fvc_reconstruct((-ghf.*fvc_snGrad(rho,xC,xF)-fvc_snGrad(p_rgh,xC,xF)).*Sf,Sf).*V;
% FRHS=fvc_reconstruct(-fvc_snGrad(p_rgh,xC,xF).*Sf,Sf).*V;

% Final assembling
UaEqnM = M + A - ASp - B + C -CSp + RHS1;  
UaEqnRHS =  RHS + ARHS - BRHS + CRHS + RHS2;

%disp('Solving for UEqna')
%Ua.internal = UaEqnM\(UaEqnRHS+FRHS); 

clear M RHS A ARHS ASp B BRHS RHS1 RHS2 C CSp;

%     {
%         volTensorField gradUbT(T(fvc::grad(Ub)));
%         volTensorField Rcb
%         (
%             "Rcb",
%             ((2.0/3.0)*I)*(k + nuEffb*tr(gradUbT)) - nuEffb*gradUbT
%         );
% 
%         surfaceScalarField phiRb
%         (
%             -fvc::interpolate(nuEffb)*mesh.magSf()*fvc::snGrad(beta)
%             /fvc::interpolate(beta + scalar(0.001))
%         );
%
phiRb = -nuEffb.*abs(Sf).*fvc_SnGrad(beta, xC, xF)./fvc_interpolate(assign(beta,constField(0.001,N),'+'), w, xC, xF);

disp('Assembling UEqnb')
%         UbEqn =
%         (
%             (scalar(1) + Cvm*rhob*alpha/rhob)*
%             (
%                 fvm::ddt(Ub)
%               + fvm::div(phib, Ub, "div(phib,Ub)")
%               - fvm::Sp(fvc::div(phib), Ub)
%             )
% 
%           - fvm::laplacian(nuEffb, Ub)
%           + fvc::div(Rcb)
% 
%           + fvm::div(phiRb, Ub, "div(phib,Ub)")
%           - fvm::Sp(fvc::div(phiRb), Ub)
% 
%           + (fvc::grad(beta)/(fvc::average(beta) + scalar(0.001)) & Rcb)
%          ==
%         //  g                          // Buoyancy term transfered to p-equation
%           - fvm::Sp(alpha/rhob*K, Ub)
%         //+ alpha/rhob*K*Ua            // Explicit drag transfered to p-equation
%           + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
%         );
% 
%         UbEqn.relax();
%     }
% }

[M, RHS] = fvm_ddt(constField(1,N),constField(1,N),Ub0,V,dt,1); % fvm::ddt(Ub)
[A, ARHS]= fvm_div_flux_cell(phib0,Ub0,xC,xF,w,1); % fvm::div(phib, Ub, "div(phib,Ub)")
ASp = fvm_Sp(fvc_div_face(phib0, V),V); % fvm::Sp(fvc::div(phib), Ub)

% (scalar(1) + Cvm*rhob*alpha/rhob)
M = diag(ones(N,1)+(transportProperties.Cvm)*(transportProperties.phasea.rho)*(alpha0.internal)/(transportProperties.phaseb.rho))*M;
A = diag(ones(N,1)+(transportProperties.Cvm)*(transportProperties.phasea.rho)*(alpha0.internal)/(transportProperties.phaseb.rho))*A;
ASp = diag(ones(N,1)+(transportProperties.Cvm)*(transportProperties.phasea.rho)*(alpha0.internal)/(transportProperties.phaseb.rho))*ASp;
RHS = (ones(N,1)+(transportProperties.Cvm)*(transportProperties.phasea.rho)*(alpha0.internal)/(transportProperties.phaseb.rho)).*RHS;
ARHS = (ones(N,1)+(transportProperties.Cvm)*(transportProperties.phasea.rho)*(alpha0.internal)/(transportProperties.phaseb.rho)).*ARHS;

% - fvm::laplacian(nuEffb, Ub)
[B, BRHS] = fvm_laplacian(Ub0,nuEffb,xC,xF,Sf);

[C, CRHS]= fvm_div_flux_cell(phiRb,Ub0,xC,xF,w,1); %  fvm::div(phiRb, Ub, "div(phib,Ub)")
CSp = fvm_Sp(fvc_div_face(phiRb, V),V);           %  fvm::Sp(fvc::div(phiRb), Ub)

%-------------- = -----------------------
% - fvm::Sp(beta/rhoa*dragCoef, Ua)
RHS1 = fvm_Sp(alpha0.internal/(transportProperties.phaseb.rho).*K, V);
RHS2 = alpha0.internal/(transportProperties.phasea.rho).*(liftCoef + (transportProperties.Cvm)*(transportProperties.phasea.rho)*DDtUa);

% Pressure terms are calculated apart in order to correctly calculate the H operator
FRHS=fvc_reconstruct((-ghf.*fvc_snGrad(rho,xC,xF)- fvc_snGrad(p_rgh,xC,xF)).*Sf,Sf).*V;

% Final assembling
UbEqnM =  M + A - ASp - B + C -CSp + RHS1;
UbEqnRHS =  RHS + ARHS - BRHS + CRHS + RHS2;

%disp('Solving for UEqnb')
%Ub.internal = UbEqnM\(UbEqnRHS+FRHS);

clear M RHS A ARHS ASp B BRHS RHS1 RHS2 FRHS C CSp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%