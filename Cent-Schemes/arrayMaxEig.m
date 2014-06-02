function [AME]=arrayMaxEig(A)
  % Gives the max eigenvalue for an array of matrices
  %
  % [AME]=arrayMaxEig(A)
  %
  % A: array of matrices
  % AME: array of max eigenvalues

  % Gets the number of matrices in array
  N=size(A,3);

  % Memory allocation
  AME=zeros(N,1);
  
  for i=1:N 
    AME(i,1)=max(max(eig(A(:,:,i))));
  end
end