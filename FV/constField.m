function [field]=constField(val,N)
  % Gives a field with full const values in internal field and BC's
  %
  % [field]=constField(value,N)
  %
  % field: the one field
  % val: const value
  % N: number of cells
  
  field.internal=ones(N,1)*val;
  
  field.left.type='V';
  field.left.value=val;
  field.left.setvalue=val;
  field.right.type='V';
  field.right.value=val;
  field.right.setvalue=val;

end