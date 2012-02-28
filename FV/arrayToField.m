function [field]=arrayToField(array)
  % Converts an array in a field with V BC given
  % by nearest cells
  %
  % [field]=arrayToField(array)
  %
  % field: the generated field
  % array: array of values
  
  field.internal=array;
  
  field.left.type='V';
  field.left.value=array(1);
  field.left.setvalue=array(1);
  field.right.type='V';
  field.right.value=array(end);
  field.right.setvalue=array(end);

end