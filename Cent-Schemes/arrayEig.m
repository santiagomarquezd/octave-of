function [arrayV,arrayVT,arrayLambda]=arrayEig(A)
  % Gives the eigenvalues matrices and eigenvectors matrices
  % for an array of matrices
  %
  % [arrayV,arrayVT,arrayLambda]=arrayEig(A)
  %
  % A: array of matrices
  % arrayV: array of eigenvector matrices
  % arrayVT: array of transposed eigenvector matrices
  % arrayLambda: array of eigenvalues matrices


  % Gets the number of matrices in array
  N=size(A,3);

  % Memory allocation
  arrayV=A;
  arrayLambda=A; 

  for i=1:N 
    [arrayV(:,:,i),arrayLambda(:,:,i)]=eig(A(:,:,i));
    arrayVT(:,:,i)=arrayV(:,:,i).';
  end
end
