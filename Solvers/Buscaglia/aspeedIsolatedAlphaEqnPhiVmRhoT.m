function [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnPhiVmRhoT(alphag,phic,rhol,rhog,V0,rhomT)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem for isolated alpha equation
    %
    % d alphag/dt+d(alphag*U+alphag*(1-cp)*Vr)/dx=0
    %
    % The eigenvalues are calculated without having rhom and Vr
    % wich can be calculated from theoretical formulas
    % rhom=rhog*alphag+(1-alphag)*rhol
    % Vr=V0*(1-alphag)
    %
    % [a_j_minus_half,a_j_plus_half]=aaspeedIsolatedAlphaEqn(alphag,rhol,rhog,V0,U)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % alphag: gas volume fraction
    % rhol: liquid density
    % rhog: gas density
    % V0: velocity constant for relative velocity law
    % U: mixture velocity

    % Notation simplification
    alphag=alphag.internal;
    rhomT=rhomT.internal;
        
    rhom=(1-alphag).*rhol+alphag.*rhog;
    Vpq=V0*(1-alphag);

    %keyboard; pause;
    eig1R=(-V0.*(1-rhog.*alphag./rhomT).*rhomT-V0.*(1-alphag).*rhomT.*rhog./rhomT).*rhog.*alphag./rhom+(phic(1:end-1)+V0.*(1-alphag).*(1-rhog./rhomT.*alphag).*rhomT).*(rhog.*rhom+rhog.*alphag.*(rhog-rhol))./(rhom.^2);
    eig1L=(-V0.*(1-rhog.*alphag./rhomT).*rhomT-V0.*(1-alphag).*rhomT.*rhog./rhomT).*rhog.*alphag./rhom+(phic(2:end)+V0.*(1-alphag).*(1-rhog./rhomT.*alphag).*rhomT).*(rhog.*rhom+rhog.*alphag.*(rhog-rhol))./(rhom.^2);
      
    aLeft=eig1R(1);
    aRight=eig1L(end);

    a_j_minus_half=[aLeft;max(eig1R(2:end),eig1L(1:end-1))];
    a_j_plus_half=[max(eig1R(2:end),eig1L(1:end-1));aRight];
	  
end	  