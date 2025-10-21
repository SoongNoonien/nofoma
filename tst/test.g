#LoadPackage("nofoma");
#SetInfoLevel(Infonofoma,1);

## Some test matrices; see also
##  A. Steel,  A new algorithm for the computation of  canonical forms of
##               matrices over fields, J. Symb. Comp. 24 (1997), 409-432.
## 

A:=[[0,0,2,2,0,2],[2,4,2,6,0,4],[3,6,0,6,1,4]];;

bev:=TransposedMat(
[[2,0,0,0,0,0,0], [2,4,1,-1,-7,-2,-1], [0,0,1,0,0,0,0], [1,0,0,1,0,0,0],
[0,0,0,0,1,0,0], [2,1,1,-1,-5,1,-1], [1,0,1,0,0,0,1]]);

steel:=
[[-23,19,-9,-75,34,9,56,15,-34,-9], [-2,2,-1,-6,3,1,4,2,-3,0], 
[4,-4,3,10,-5,-1,-6,-4,5,1], [-2,2,-1,-5,3,1,3,2,-3,0], 
[0,0,0,0,2,0,0,0,0,0], [12,-12,6,33,-18,-4,-18,-12,18,0], 
[-1,-3,0,2,1,0,1,1,2,1], [-26,22,-10,-83,36,10,61,18,-39,-10], 
[-1,-3,0,1,1,0,2,1,2,0], [8,-12,4,27,-12,-4,-12,-7,15,0]];;

low:=Z(19)^0 *[[16,0,0,0,0,0,0,0,0,0],[5,2,0,0,0,0,0,0,0,0],
[9,6,7,0,0,0,0,0,0,0],[15,9,1,2,0,0,0,0,0,0],[7,14,8,2,9,0,0,0,0,0],
[12,3,17,5,1,12,0,0,0,0],[7,11,0,4,6,9,10,0,0,0],[8,3,15,16,17,18,18,12,0,0],
[6,3,7,12,1,12,11,14,10,0],[18,14,7,17,16,15,13,13,3,8]];;

#a:=(RandomLowerTriangularMat(1000,1000,GF(4)));;
#nma:=NeuMaximalVectorMat(a);;time;
#Degree(nma[2])=Length(SpinMatVector(a,nma[1])[2]);
#nma[2]=MinimalPolynomial(a);

nfmJordanBlock:=function(n,c)
  local a,i;
  a:=IdentityMat(n,DefaultField(c));
  for i in [1..n-1] do
    a[i][i]:=c;
    a[i+1][i]:=c^0;
  od;
  a[n][n]:=c;
  return a;
end;

NPmatrixM3:=function(n)
  local a,b,c,i;
  if n=600 then 
    a:=BuildBlockDiagonalMat(nfmJordanBlock(300,Z(5)),
                             Z(5)*IdentityMat(300,GF(5)));
    b:=Random(GL(600,GF(5)));
    return b*a*b^-1;
  elif n=1200 then 
    a:=nfmJordanBlock(2,Z(3));
    c:=nfmJordanBlock(2,Z(3));
    for i in [2..400] do 
      c:=BuildBlockDiagonalMat(c,a);
    od;
    c:=BuildBlockDiagonalMat(c,Z(3)*IdentityMat(400,GF(3)));
    b:=Random(GL(1200,GF(3)));
    return b*c*b^-1;
  fi;
end;

mat1:=function(mat)
  local a,a1,b,i;
  a1:=TransposedMat(Concatenation(mat,mat));
  a:=[];
  for i in [1..Length(a1)-1] do
    Add(a,a1[i]);
  od;
  Add(a,0*a1[1]);
  b:=TransposedMat(Concatenation(TransposedMat(mat),TransposedMat(mat)));
  return Concatenation(a,b);
end;

diagmat:=function(l)
  local m,i;
  m:=l[1]*IdentityMat(Length(l));
  for i in [2..Length(l)] do
    m[i][i]:=l[i];
  od;
 return m;
end;

ddd:=diagmat([1,2,-1,0,4,3,3,4,5,6,7,-1,5,4,0,0,3,2,1]);

# Tests
test:=function(a1)
  local a,aa;
  a:=CyclicChainMat(a1);
  a:=MinPolyMat(a1);
  a:=InvariantFactorsMat(a1);
  a:=FrobeniusNormalForm(a1);
  if CheckFrobForm(a1,a)=false
    then Error(" --> tests not OK !\n");
  fi;
  if CheckJordanChev(a1,JordanChevalleyDecMat(a1,MinPolyMat(a1)))=false
    then Error(" --> tests not OK !\n");
  fi;
end;

FrobTest:=function(lev)
  local f;
  test(steel); test(bev); test(ddd);
  test(Z(4)^0*steel); test(Z(4)^0*bev); test(Z(4)^0*ddd);
  test(mat1(steel)); test(Z(29)*mat1(steel));
  test(RandomMat(15,15,Integers));
  test(RandomMat(20,20,GF(49)));
  #test(RandomLowerTriangularMat(10,10,GF(19)));
  Print("\n--> tests OK !\n\n");
end;

