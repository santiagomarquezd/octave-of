function [AMN]=arrayMatrixNorm(A,n)
  % Gives the norm for an array of matrices
  %
  % [AMN]=arrayMatrixNorm(A)
  %
  % A: array of matrices
  % n: type of norm (1, Inf, 2, etc)
  % AMN: array of norms

  % Gets the number of matrices in array
  N=size(A,3);

  % Memory allocation
  AMN=zeros(N,1);
  
  for i=1:N 
    AMN(i,1)=norm(A(:,:,i),n);
  end
end