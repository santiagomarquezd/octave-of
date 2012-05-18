function [AInv]=matrixArrayInverse(A)
  % Gives the inverse of an array of square matrices
  %
  % [AInv]=matrixArrayInverse(A)
  %
  % AInv: array of matrix inverses
  % A: array of matrices


  % Gets the number of matrices in the arrays
  N=size(A,3);
  
  % Memory allocation
  AInv=A; 
  
  for i=1:N 
    AInv(:,:,i)=inv(A(:,:,i));
  end
end