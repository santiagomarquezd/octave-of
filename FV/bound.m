function [out]=bound(in,op,val)
  % Bounds an array by op and val
  %
  % [out]=bound(out,op,val)
  %
  % out: resulting array
  % in: entering array
  % op: operation ('max' or 'min')
  % val: bounding value
  
  % General assignment
  out=in;
  
  % Bounding
  if (op=='max')
    out.internal(in.internal>val)=val;  
  elseif (op=='min')
    out.internal(in.internal<val)=val;  
  end
  
end