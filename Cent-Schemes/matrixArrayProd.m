function [AprodB]=matrixArrayProd(A,B)
  % Gives the matrix product of two array of matrices
  %
  % [AprodB]=matrixArrayProd(A,B)
  %
  % AprodB: array of matrix products
  % A: array of first factors
  % B: array of second factors

  % Gets the number of matrices in the arrays
  NA=size(A,3);
  NB=size(B,3);
  
  % Size checking
  if (NA!=NB) 
    disp('Non conformant size of arrays');
    return;
  elseif (size(A(:,:,1),1)!=size(B(:,:,1),1)||size(A(:,:,1),2)!=size(B(:,:,1),2))
    disp('Non conformant size of matrices');
    return;
  end

  % Memory allocation
  AprodB=A; 
  
  for i=1:NA 
    AprodB(:,:,i)=A(:,:,i)*B(:,:,i);
  end
end