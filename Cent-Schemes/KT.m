function [u,v]=KT(u0,v0,uLimited,vLimited,ufluxfunction,vfluxfunction,Um,a,dx,dt,V0,rhol,rhog,aexp)
    % Advances u and v fields in time by means of Kurganov & Tadmor scheme
    %
    % [u,v]=KT(u0,v0,uLimited,vLimited,ufluxfunction,vfluxfunction,Um,a,dx,dt,V0,rhol,rhog,aexp)
    %
    % u0: actual value of u
    % v0: actual value of v
    % u: first field advanced in time
    % v: second field advanced in time
    % uLimited: limited values for u at actual tiem
    % vLimited: limited values for u at actual tiem
    % ufluxfunction: flux function for u
    % vfluxfunction: flux function for v
    % Um: linear advection field
    % a: local velocities at interfaces (eigenvalues)
    % dx: spatial step
    % dt: time step
    % V0: constant for flux function
    % rhol: liquid density
    % rhog: gas density
    % aexp: exponent in advective velocity calculation
    
    % Allocation
    u=u0;
    v=v0;
    
    % Pre-processing
    lambda=dt/dx;

    %keyboard; pause;
    
    % Star fluxes calculation (see http://en.wikipedia.org/wiki/MUSCL_scheme)
    Fu_star_minus_half=1/2*(ufluxfunction(uLimited.u_i_m_h_r,Um,V0,rhol,rhog,aexp)+...
						ufluxfunction(uLimited.u_i_m_h_l,Um,V0,rhol,rhog,aexp)-...
						a(1:end-1).*(uLimited.u_i_m_h_r-uLimited.u_i_m_h_l));
    Fu_star_plus_half=1/2*(ufluxfunction(uLimited.u_i_p_h_r,Um,V0,rhol,rhog,aexp)+...
						ufluxfunction(uLimited.u_i_p_h_l,Um,V0,rhol,rhog,aexp)-...
						a(2:end).*(uLimited.u_i_p_h_r-uLimited.u_i_p_h_l));
    Fv_star_minus_half=1/2*(vfluxfunction(vLimited.u_i_m_h_r,Um,V0,rhol,rhog,aexp)+...
						vfluxfunction(vLimited.u_i_m_h_l,Um,V0,rhol,rhog,aexp)-...
						a(1:end-1).*(vLimited.u_i_m_h_r-vLimited.u_i_m_h_l));
    Fv_star_plus_half=1/2*(vfluxfunction(vLimited.u_i_p_h_r,Um,V0,rhol,rhog,aexp)+...
						vfluxfunction(vLimited.u_i_p_h_l,Um,V0,rhol,rhog,aexp)-...
						a(2:end).*(vLimited.u_i_p_h_r-vLimited.u_i_p_h_l));
    % Impermeable walls BC's
    if 0
      Fu_star_minus_half(1)=0;
      Fu_star_plus_half(end)=0;
      Fv_star_minus_half(1)=0;
      Fv_star_plus_half(end)=0;
    end

    % Time advancement
    u.internal=u0.internal-lambda*(Fu_star_plus_half-Fu_star_minus_half);
    v.internal=v0.internal-lambda*(Fv_star_plus_half-Fv_star_minus_half);
    %end
end
