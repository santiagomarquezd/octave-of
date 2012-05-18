function [AAbs]=matrixArrayAbs(V,lambda)
  % Gives the absolute value of an array of square matrices
  %
  % [AAbs]=matrixArrayAbs(A,V)
  %
  % AAbs: array of matrix absolute values
  % A: array of matrices
  % V: array of eigenvectors
  % lambda: array of matrices eigenvalues


  % Inverse of eigenvector matrices
  VInv=matrixArrayInverse(V);

  % Array of matrices of absolute values of eigenvalues
  N=size(V,3);
  AbsLambda=zeros(2,2,N);
  for i=1:N
    AbsLambda(:,:,i)=diag(abs(diag(lambda(:,:,i))));
  end
 

  % Calculation
  AAbs=matrixArrayProd(matrixArrayProd(V,AbsLambda),VInv);

end