function [uLimited]=limitedValues(u,phi,dx,dt)
    % Gives the limited values for +1/2*dx and -1/2*dx interfaces
    % at both sides (L, R) needed for Kurganov & Tadmor scheme
    %
    % [ulimited]=limitedValues(u,phi,dx,dt)
    %
    % uLimited.u_i_p_h_l: u in i+1/2,L
    % uLimited.u_i_p_h_r: u in i+1/2,R
    % uLimited.u_i_m_h_l: u in i-1/2,L
    % uLimited.u_i_m_h_r: u in i-1/2,R
    % u: field values at cell centres
    % phi: Sweby's function values
    % dx: spatial step
    % dt: time-step
    %
    % **************************************************************
    %                           CAUTION
    % **************************************************************
    %
    % Fixed value BC at left and zero gradient BC at right are assumed
    % for u field
    
    % Precomputing
    lambda=dt/dx;
    
    % Memory allocation
    uLimited.u_i_m_h_r=uLimited.u_i_m_h_l=uLimited.u_i_p_h_r=uLimited.u_i_p_h_l=0*u.internal;
    
    % Limited values at interfaces
    uLimited.u_i_p_h_l(1:end-1)=u.internal(1:end-1)+0.5*phi(1:end-1).*(u.internal(2:end)-u.internal(1:end-1));
    % Zero gradient BC
    uLimited.u_i_p_h_l(end)=u.internal(end);
   
    uLimited.u_i_p_h_r(1:end-2)=u.internal(2:end-1)-0.5*phi(2:end-1).*(u.internal(3:end)-u.internal(2:end-1));
    uLimited.u_i_p_h_r(end-1)=u.internal(end)-0.5*phi(end)*(u.right.setvalue-u.internal(end));
    % Zero gradient BC
    uLimited.u_i_p_h_r(end)=u.right.setvalue;
  
    uLimited.u_i_m_h_l(2:end)=u.internal(1:end-1)+0.5*phi(1:end-1).*(u.internal(2:end)-u.internal(1:end-1));
    % Fixed value BC
    uLimited.u_i_m_h_l(1)=u.left.setvalue;
    
    uLimited.u_i_m_h_r(1:end-1)=u.internal(1:end-1)-0.5*phi(1:end-1).*(u.internal(2:end)-u.internal(1:end-1));
    % Zero gradient BC
    uLimited.u_i_m_h_r(end)=u.internal(end);
end

