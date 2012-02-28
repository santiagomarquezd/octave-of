function [u,v]=KT(u0,v0,uLimited,vLimited,ufluxfunction,vfluxfunction,a,dx,dt,Vr,rhol,rhog,model)
    % Advances u and v fields in time by means of Kurganov & Tadmor scheme
    %
    % [u,v]=KT(uLimited,vLimited,ufluxfunction,vfluxfunction)
    %
    % u0: actual value of u
    % v0: actual value of v
    % u: first field advanced in time
    % v: second field advanced in time
    % uLimited: limited values for u at actual tiem
    % vLimited: limited values for u at actual tiem
    % ufluxfunction: flux function for u
    % vfluxfunction: flux function for v
    % a: local velocities (eigenvalues)
    % dx: spatial step
    % dt: time step
    % Vr: relative velocity
    % rhol: liquid density
    % rhog: gas density
    % model: 1. No Vr; 2: with Vr
    
    % Allocation
    u=u0;
    v=v0;
    
    % Pre-processing
    lambda=dt/dx;
    
    % Star fluxes calculation (see http://en.wikipedia.org/wiki/MUSCL_scheme)
    Fu_star_minus_half=1/2*(ufluxfunction(uLimited.u_i_m_h_r,vLimited.u_i_m_h_r,Vr,rhol,rhog,model)+...
						ufluxfunction(uLimited.u_i_m_h_l,vLimited.u_i_m_h_l,Vr,rhol,rhog,model)-...
						a(1:end-1).*(uLimited.u_i_m_h_r-uLimited.u_i_m_h_l));
    Fu_star_plus_half=1/2*(ufluxfunction(uLimited.u_i_p_h_r,vLimited.u_i_p_h_r,Vr,rhol,rhog,model)+...
						ufluxfunction(uLimited.u_i_p_h_l,vLimited.u_i_p_h_l,Vr,rhol,rhog,model)-...
						a(2:end).*(uLimited.u_i_p_h_r-uLimited.u_i_p_h_l));
	Fv_star_minus_half=1/2*(vfluxfunction(uLimited.u_i_m_h_r,vLimited.u_i_m_h_r,Vr,rhol,rhog,model)+...
						vfluxfunction(uLimited.u_i_m_h_l,vLimited.u_i_m_h_l,Vr,rhol,rhog,model)-...
						a(1:end-1).*(vLimited.u_i_m_h_r-vLimited.u_i_m_h_l));
    Fv_star_plus_half=1/2*(vfluxfunction(uLimited.u_i_p_h_r,vLimited.u_i_p_h_r,Vr,rhol,rhog,model)+...
						vfluxfunction(uLimited.u_i_p_h_l,vLimited.u_i_p_h_l,Vr,rhol,rhog,model)-...
						a(2:end).*(vLimited.u_i_p_h_r-vLimited.u_i_p_h_l));
    

	% Time advancement
	u.internal=u0.internal-lambda*(Fu_star_plus_half-Fu_star_minus_half);
	v.internal=v0.internal-lambda*(Fv_star_plus_half-Fv_star_minus_half);
end
