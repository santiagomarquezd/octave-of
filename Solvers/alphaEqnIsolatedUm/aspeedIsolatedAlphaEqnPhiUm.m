function [a_j_minus_half,a_j_plus_half]=aspeedIsolatedAlphaEqnPhiUm(alphalLimited,phic,rhol,rhog,V0,a)
    % Gives the max of spectral radius at interfaces for 
    % Kurganov and Tadmor's kind of schemes
    % in simplified mixture problem for isolated alpha equation
    % The center of volume velocity, Um,  is given at faces as a flux
    % (this case is simply since the cross areas is 1)
    %
    % d alphal/dt+d(alphal*Um-alphal^a*(1-alphal)*Urlg)/dx=0
    %
    % [a_j_minus_half,a_j_plus_half]=aaspeedIsolatedAlphaEqn(alphal,phic,rhol,rhog,V0,a)
    %
    % a_j_minus_half: max of spectral radius at left interfaces
    % a_j_plus_half: max of spectral radius at right interfaces
    % alphalLimited: limited liquid volume fraction given at i+1/2,l,r and i-1/2,l,r
    % phic: face flux for the velocity of center of volume
    % rhol: liquid density
    % rhog: gas density
    % V0: velocity constant for relative velocity law
    % a: exponent constant for relative velocity law

    % uLimited structure
    % uLimited.u_i_m_h_l
    % uLimited.u_i_m_h_r
    % uLimited.u_i_p_h_l
    % uLimited.u_i_p_h_r

    % Eigenvalues calculation
    eig_i_m_h_l=-V0.*((a+1).*(alphalLimited.u_i_m_h_l).^a-(a+2).*alphalLimited.u_i_m_h_l.^(a+1))+phic(1:end-1);
    eig_i_m_h_r=-V0.*((a+1).*(alphalLimited.u_i_m_h_r).^a-(a+2).*alphalLimited.u_i_m_h_r.^(a+1))+phic(1:end-1);
    eig_i_p_h_l=-V0.*((a+1).*(alphalLimited.u_i_p_h_l).^a-(a+2).*alphalLimited.u_i_p_h_l.^(a+1))+phic(2:end);
    eig_i_p_h_r=-V0.*((a+1).*(alphalLimited.u_i_p_h_r).^a-(a+2).*alphalLimited.u_i_p_h_r.^(a+1))+phic(2:end);

    % 
    a_j_minus_half=max(abs(eig_i_m_h_l),abs(eig_i_m_h_r));
    a_j_plus_half=max(abs(eig_i_p_h_l),abs(eig_i_p_h_r));
	  
end	  



