function [char1,char2]=primToChar(prim1,prim2,eigens,rhol,rhog)
  % Transforms primitive variables to characteristic ones
  % via inverse eigenvectors matrix
  % (Characteristic velocities MUST retain the direction
  % of original velocities for correct BC setting)
  %
  % [char1,char2]=primToChar(prim1,prim2,eigens)
  % 
  % char1,char2: output has characteristic variables
  % prim1,prim2: primitive variables
  % eigens: function to find the eigenvectors matrix
  % rhol: liquid density
  % rhog: gas density
  
  % Variable's size
  N=size(prim1.internal,1);
  
  % Memory allocation
  char1=prim1;
  char2=prim2;
  
  % Transformation 
  % Internal field
  for i=1:size
    VV=eigens(prim1.internal(i),prim2.internal(i),rhol,rhog);
    VVinv=inverse(VV);
    c=VVinv*[prim1.internal(i);prim2.internal(i)];
    char1=c(1,1);
    char2=c(2,1);  
  end
  
  % BC's
  if (prim1.left.type=='V' && prim2.left.type=='V')
    VV=eigens(prim1.left.setvalue(i),prim2.left.setvalue(i),rhol,rhog);
    VVinv=inverse(VV);
    c=VVinv*[prim1.left.setvalue(i);prim2.left.setvalue(i)];
    char1.left.value=c(1,1);
    char2.left.value=c(2,1);
  elseif (prim1.right.type=='V' && prim2.right.type=='V')
    VV=eigens(prim1.right.setvalue(i),prim2.right.setvalue(i),rhol,rhog);
    VVinv=inverse(VV);
    c=VVinv*[prim1.right.setvalue(i);prim2.right.setvalue(i)];
    char1.right.value=c(1,1);
    char2.right.value=c(2,1);
  end
  

end
