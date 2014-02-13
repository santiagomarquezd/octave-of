function [r]=rvalue(u,estab)
    % r variable calculation for TVD limiter
    %
    % [r]=rvalue(u)
    %
    % r: r field
    % u: given field
    % estab: estabilization factor for quotients
    
    % Memory allocation
    r=u.internal*0;
    
    % Internal
    r(2:end-1)=(u.internal(2:end-1)-u.internal(1:end-2))./(u.internal(3:end)-u.internal(2:end-1)+estab);
    % BC's
    r(1)=2*(u.internal(1)-u.left.setvalue)/(u.internal(2)-u.internal(1)+estab);
    r(end)=(u.internal(end)-u.internal(end-1))/(u.right.setvalue-u.internal(end)+estab)/2;
    
end



 
  
