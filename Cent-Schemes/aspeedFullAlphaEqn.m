function [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(alphag,rhol,rhog,V0)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem
    %
    % dU/dt+d(rhom*U*U)/dx=S(dp/dx-rhom*g) (negligible drift term)
    % d alphag/dt+d(alphag*U+(1-cp)*Vr)/dx=0
    %
    % The eigenvalues are calculated without having rhom and U
    % wich can be calculated from theoretical formulas
    % rhom=rhog*alphag+(1-alphag)*rhol
    % Vr=V0*(1-alphag)
    % U=(rhog/rhom-1)*alphag*Vr
    %
    % [a_j_minus_half,a_j_plus_half]=aspeedFullAlphaEqn(U,alphag,rhol,rhog,V0)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % U: mixture velocity
    % rhol: liquid density
    % rhog: gas density
    % alphag: gas volume fraction
    
	% Alpha is flatten for an easier calculus (BC's and internal at the same time)
	alphag=[alphag.left.setvalue; alphag.internal; alphag.right.setvalue];
	
	eig1=-(sqrt((1-alphag).^2.*alphag.^2.*(1000000.*(15.*alphag.^6-20.*alphag.^5+6.*alphag.^4)+1000000000000.*(15.*alphag.^6-40.*alphag.^5+36.*alphag.^4-12.*alphag.^3+alphag.^2)+1000000000000000000.*(alphag.^6-4.*alphag.^5+6.*alphag.^4-4.*alphag.^3+alphag.^2)+alphag.^6+1000000000000000.*(-6.*alphag.^6+20.*alphag.^5-24.*alphag.^4+12.*alphag.^3-2.*alphag.^2)+1000.*(4.*alphag.^5-6.*alphag.^6)+1000000000.*(-20.*alphag.^6+40.*alphag.^5-24.*alphag.^4+4.*alphag.^3)+(alphag+1000.*(1-alphag)).*(1000000000.*(20.*alphag.^5-48.*alphag.^4+36.*alphag.^3-8.*alphag.^2)+1000.*(10.*alphag.^5-8.*alphag.^4)+1000000000000000.*(2.*alphag.^5-8.*alphag.^4+12.*alphag.^3-8.*alphag.^2+2.*alphag)-2.*alphag.^5+1000000000000.*(-10.*alphag.^5+32.*alphag.^4-36.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000.*(-20.*alphag.^5+32.*alphag.^4-12.*alphag.^3))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(6.*alphag.^4-12.*alphag.^3+6.*alphag.^2)+1000000000000.*(alphag.^4-4.*alphag.^3+6.*alphag.^2-4.*alphag+1)+alphag.^4+1000000000.*(-4.*alphag.^4+12.*alphag.^3-12.*alphag.^2+4.*alphag)+1000.*(4.*alphag.^3-4.*alphag.^4))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).^2.*V0.^2+0.282.*(1-alphag).*alphag.*((alphag+1000.*(1-alphag)).*(1000000.*(16.*alphag.^6-30.*alphag.^5+16.*alphag.^4-2.*alphag.^3)+1000000000000.*(16.*alphag.^6-58.*alphag.^5+80.*alphag.^4-52.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000000000000.*(-4.*alphag.^6+18.*alphag.^5-32.*alphag.^4+28.*alphag.^3-12.*alphag.^2+2.*alphag)+1000.*(4.*alphag.^5-4.*alphag.^6)+1000000000.*(-24.*alphag.^6+66.*alphag.^5-64.*alphag.^4+26.*alphag.^3-4.*alphag.^2))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(12.*alphag.^5-26.*alphag.^4+16.*alphag.^3-2.*alphag.^2)+1000000000000.*(4.*alphag.^5-18.*alphag.^4+32.*alphag.^3-28.*alphag.^2+12.*alphag-2)+1000.*(4.*alphag.^4-4.*alphag.^5)+1000000000.*(-12.*alphag.^5+40.*alphag.^4-48.*alphag.^3+24.*alphag.^2-4.*alphag))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.079524.*(alphag+1000.*(1-alphag)).^2.*(1000000.*(4.*alphag.^6-8.*alphag.^5+4.*alphag.^4)+1000000000000.*(4.*alphag.^6-20.*alphag.^5+41.*alphag.^4-44.*alphag.^3+26.*alphag.^2-8.*alphag+1)+1000000000.*(-8.*alphag.^6+28.*alphag.^5-36.*alphag.^4+20.*alphag.^3-4.*alphag.^2)))+(1-alphag).*alphag.*(1000000.*(3.*alphag.^3-4.*alphag.^2+alphag)+alphag.^3+1000000000.*(-alphag.^3+2.*alphag.^2-alphag)+1000.*(2.*alphag.^2-3.*alphag.^3)+(alphag+1000.*(1-alphag)).*(1000.*(6.*alphag.^2-6.*alphag)-3.*alphag.^2+1000000.*(-3.*alphag.^2+6.*alphag-3))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.282.*(alphag+1000.*(1-alphag)).*(1000000.*(2.*alphag.^3-5.*alphag.^2+4.*alphag-1)+1000.*(2.*alphag.^2-2.*alphag.^3)))./((alphag+1000.*(1-alphag)).*(1000000.*(2.*alphag.^2-4.*alphag+2)+2.*alphag.^2+1000.*(4.*alphag-4.*alphag.^2)));
	
	eig2=(sqrt((1-alphag).^2.*alphag.^2.*(1000000.*(15.*alphag.^6-20.*alphag.^5+6.*alphag.^4)+1000000000000.*(15.*alphag.^6-40.*alphag.^5+36.*alphag.^4-12.*alphag.^3+alphag.^2)+1000000000000000000.*(alphag.^6-4.*alphag.^5+6.*alphag.^4-4.*alphag.^3+alphag.^2)+alphag.^6+1000000000000000.*(-6.*alphag.^6+20.*alphag.^5-24.*alphag.^4+12.*alphag.^3-2.*alphag.^2)+1000.*(4.*alphag.^5-6.*alphag.^6)+1000000000.*(-20.*alphag.^6+40.*alphag.^5-24.*alphag.^4+4.*alphag.^3)+(alphag+1000.*(1-alphag)).*(1000000000.*(20.*alphag.^5-48.*alphag.^4+36.*alphag.^3-8.*alphag.^2)+1000.*(10.*alphag.^5-8.*alphag.^4)+1000000000000000.*(2.*alphag.^5-8.*alphag.^4+12.*alphag.^3-8.*alphag.^2+2.*alphag)-2.*alphag.^5+1000000000000.*(-10.*alphag.^5+32.*alphag.^4-36.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000.*(-20.*alphag.^5+32.*alphag.^4-12.*alphag.^3))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(6.*alphag.^4-12.*alphag.^3+6.*alphag.^2)+1000000000000.*(alphag.^4-4.*alphag.^3+6.*alphag.^2-4.*alphag+1)+alphag.^4+1000000000.*(-4.*alphag.^4+12.*alphag.^3-12.*alphag.^2+4.*alphag)+1000.*(4.*alphag.^3-4.*alphag.^4))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).^2.*V0.^2+0.282.*(1-alphag).*alphag.*((alphag+1000.*(1-alphag)).*(1000000.*(16.*alphag.^6-30.*alphag.^5+16.*alphag.^4-2.*alphag.^3)+1000000000000.*(16.*alphag.^6-58.*alphag.^5+80.*alphag.^4-52.*alphag.^3+16.*alphag.^2-2.*alphag)+1000000000000000.*(-4.*alphag.^6+18.*alphag.^5-32.*alphag.^4+28.*alphag.^3-12.*alphag.^2+2.*alphag)+1000.*(4.*alphag.^5-4.*alphag.^6)+1000000000.*(-24.*alphag.^6+66.*alphag.^5-64.*alphag.^4+26.*alphag.^3-4.*alphag.^2))+(alphag+1000.*(1-alphag)).^2.*(1000000.*(12.*alphag.^5-26.*alphag.^4+16.*alphag.^3-2.*alphag.^2)+1000000000000.*(4.*alphag.^5-18.*alphag.^4+32.*alphag.^3-28.*alphag.^2+12.*alphag-2)+1000.*(4.*alphag.^4-4.*alphag.^5)+1000000000.*(-12.*alphag.^5+40.*alphag.^4-48.*alphag.^3+24.*alphag.^2-4.*alphag))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.079524.*(alphag+1000.*(1-alphag)).^2.*(1000000.*(4.*alphag.^6-8.*alphag.^5+4.*alphag.^4)+1000000000000.*(4.*alphag.^6-20.*alphag.^5+41.*alphag.^4-44.*alphag.^3+26.*alphag.^2-8.*alphag+1)+1000000000.*(-8.*alphag.^6+28.*alphag.^5-36.*alphag.^4+20.*alphag.^3-4.*alphag.^2)))+(1-alphag).*alphag.*(1000.*(3.*alphag.^3-2.*alphag.^2)+1000000000.*(alphag.^3-2.*alphag.^2+alphag)-alphag.^3+1000000.*(-3.*alphag.^3+4.*alphag.^2-alphag)+(alphag+1000.*(1-alphag)).*(1000000.*(3.*alphag.^2-6.*alphag+3)+3.*alphag.^2+1000.*(6.*alphag-6.*alphag.^2))).*(rhog./((1-alphag).*rhol+alphag.*rhog)-1).*V0+0.282.*(alphag+1000.*(1-alphag)).*(1000.*(2.*alphag.^3-2.*alphag.^2)+1000000.*(-2.*alphag.^3+5.*alphag.^2-4.*alphag+1)))./((alphag+1000.*(1-alphag)).*(1000000.*(2.*alphag.^2-4.*alphag+2)+2.*alphag.^2+1000.*(4.*alphag-4.*alphag.^2)));
	  
	% Local velocities at internal field
	if 0
	  a=max(abs(eig1),abs(eig2));
	else
	   a=max(abs(eig2),abs(eig2));
	end
	  
	aLeft=a(1);
	aRight=a(end);
	
	% a is turned back to internal values only
	a=a(2:end-1);
    
	a_j_minus_half=[aLeft;max(a(1:end-1),a(2:end))];
	a_j_plus_half=[max(a(1:end-1),a(2:end));aRight];
	  
end	  