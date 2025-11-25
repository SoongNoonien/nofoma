#########################################################################
##
#A  nofoma.gd                                                Meinolf Geck
##
#Y  Copyright (C) 2019/22  Lehrstuhl fuer Algebra, Universitaet Stuttgart
##
##  This is a package for computing maximal vectors, minimal polynomials, 
##  the rational canonical form (or Frobenius normal form)  and also  the
##  Jordan-Chevalley decomposition  of a  matrix  over any  field that is 
##  available in  GAP. The algorithms are based on, and a combination of:  
## 
##  K. Bongartz,  A direct approach to the rational normal form, preprint
##                                          available at arXiv:1410.1683.
## 
##  M. Neunhoeffer and  C. E. Praeger,  Computing minimal polynomials  of
##                     matrices, LMS J. Comput. Math. 11 (2008), 252-279;
## 
##  M. Geck,  On Jacob's construction  of  the  rational  canonical form, 
##                        Electron. J. Linear Algebra 36 (2020), 177-182.
##  
##  M. Geck, On the Jordan-Chevalley decomposition of a matrix, preprint.
##
##  Stuttgart, June 1, 2022.
#########################################################################

#! @Chapter The nofoma package
#! @ChapterLabel The nofoma package
#! @Section Maximal vectors and normal forms
#! Let <M>K</M> be a field and <M>A</M> be an <M>n\times n</M>-matrix over 
#! <M>K</M>. This package provides functions for computing a maximal vector
#! for <M>A</M>, the Frobenius normal form of <M>A</M> and the 
#! Jordan-Chevalley decomposition of <M>A</M>.
#! 
#! For any (row) vector <M>v\in K^n</M>, the local minimal polynomial is the 
#! unique monic polynomial <M>f\in K[X]</M> of smallest possible degree such 
#! that <M>v.f(A)=0</M>. (Note that, as usual in GAP, matrices act on the 
#! right on row vectors.) It is known that there always exists a 
#! <M>v\in K^n</M> such that the local minimal polynomial of <M>v</M> equals 
#! the minimal polynomial of <M>A</M>; such a <M>v</M> is called a 'maximal 
#! vector' for <M>A</M>. 
#! 
#! The currently best algorithm for computing the minimal polynomial of 
#! <M>A</M> is probably that of Neunhoeffer--Praeger (already implemented 
#! in GAP). One of the purposes of this package is to modify that algorithm 
#! so that it also includes the computation of a maximal vector for <M>A</M>; 
#! this is done by the function 'MaximalVectorMat'. Once this is available, 
#! the Frobenius normal form (or rational canonical form) of <M>A</M> can be 
#! computed efficiently by a recursive algorithm. 
#! 
#! Finally, we also provide a function for computing the Jordan-Chevalley
#! decomposition of <M>A</M>. Usually, this is obtained as a consequence of 
#! the Jordan normal form of <M>A</M>, but in general the latter is difficult 
#! to compute because one needs to know the eigenvalues of <M>A</M>. However,
#! there is an elegant algorithmic approach (going back to Chevalley), which 
#! is inspired by the Newton iteration for finding the zeroes of a function;
#! and it does not require the knowledge of the eigenvalues of <M>A</M>. This
#! is implemented in the function 'JordanChevalleyDecMat'.
#! 
#! Further references are provided within the help menu of each function.
#! 
#! Any comments welcome!                            Meinolf Geck, June 2022
#! 
#! @Section Copyright and installation of the nofoma package
#! The nofoma package is free software; you can redistribute it and/or 
#! modify it under the terms of the GNU General Public License as published 
#! by the Free Software Foundation; either version 2 of the License, or (at 
#! your option) any later version.
#!
#! To install this package first unpack it inside some GAP root directory
#! in the subdirectory `pkg` (see the section 'Installing a GAP Package' of
#! the GAP reference manual). Then 'nofoma' can already be loaded and used
#! (just type `LoadPackage("nofoma");`).
#!
#! @Section The main functions

DeclareInfoClass("Infonofoma"); 
##  Since 'Info' does unwanted formatting, we use (from EDIM package):
DeclareGlobalFunction("IInofoma");
DeclareGlobalFunction("nfmCoeffsPol"); 
DeclareGlobalFunction("nfmPolCoeffs");
DeclareGlobalFunction("nfmGcd");
DeclareGlobalFunction("nfmLcm");

#! @Arguments a,b
#! @Description 
#! Computes a divisor <M>a_1</M> of the polynomial <A>a</A> and a
#! divisor <M>b_1</M> of the polynomial <A>b</A> such that <M>a_1</M> and <M>b_1</M> are coprime
#! and the lcm of <A>a</A>, <A>b</A> is <M>a_1</M>*<M>b_1</M>.  This is based on Lemma 5 in <Cite Key ="Bon14"/>.
#! (see also Lemma 4.3 in <Cite Key ="Gec20"/>).
#!
#! (Note that it does not use the prime factorisation of polynomials but 
#! only gcd computations.)
#! 
#! @BeginExampleSession
#! gap> a:=x^2*(x-1)^3*(x-2)*(x-3);
#! x^7-8*x^6+24*x^5-34*x^4+23*x^3-6*x^2
#! gap> b:=x^2*(x-1)^2*(x-2)^4*(x-4);
#! x^9-14*x^8+81*x^7-252*x^6+456*x^5-480*x^4+272*x^3-64*x^2
#! gap> GcdCoprimeSplit(a,b);
#! [ x^5-4*x^4+5*x^3-2*x^2,              
#! x^4-6*x^3+12*x^2-10*x+3,               
#! x^7-12*x^6+56*x^5-128*x^4+144*x^3-64*x^2 ]  # b1
#! @EndExampleSession
DeclareGlobalFunction("GcdCoprimeSplit");

#! @Arguments A,pol,v
#! @Description
#! Returns  the row vector  obtained  by multiplying 
#! the row vector <A>v</A> with the matrix <A>pol</A>(<A>A</A>), where <A>pol</A> is the list 
#! of coefficients of a polynomial.
#!
#! @BeginExampleSession
#! gap> A:=[ [ 0, 1, 0, 1 ],
#! gap>       [ 0, 0, 0, 0 ],
#! gap>       [ 0, 1, 0, 1 ],
#! gap>       [ 1, 1, 1, 1 ] ];;
#! gap> f:=x^6-6*x^5+12*x^4-10*x^3+3*x^2;;
#! gap> v:=[ 1, 1, 1, 1];;
#! gap> l:=nfmCoeffsPol(f);
#! gap> [ 0, 0, 3, -10, 12, -6, 1 ]
#! gap> PolynomialToMat(A,last,v);
#! [ 8, -16, 8, -16 ]
#! @EndExampleSession
DeclareGlobalFunction("PolynomialToMatVec");

DeclareGlobalFunction("PolynomialToMat");

#! @Arguments A,v1,v2,pol1,pol2
#! @Description
#!  Returns,  given  a matrix  <A>A</A>,  vectors <A>v1</A>,
#!  <A>v2</A> with minimal polynomials <A>pol1</A>, <A>pol2</A>,  a new pair [<M>v</M>,<M>pol</M>],  
#!  where <M>v</M> has minimal polynomial <M>pol</M>, and <M>pol</M> is the least common
#!  multiple of <A>pol1</A> and <A>pol2</A>.  
#!  This crucially relies on  'GcdCoprimeSplit' to avoid  factorisation of 
#!  polynomials.
DeclareGlobalFunction("LcmMaximalVectorMat");

DeclareGlobalFunction("SpinMatVector1");

#! @Arguments A,v
#! @Description
#!  Computes  the  smallest subspace containing the vector  
#!  <A>v</A> that is invariant under the matrix <A>A</A>. The  output is a 
#!  quadruple, with 
#!  * 1st component = basis of that subspace in row echelon form; 
#!  * 2nd component = matrix  with  rows <A>v</A>, <A>v.A</A>, <A>v.A^2</A>, 
#!     ..., <A>v.A^{{d-1}}</A> (where <M>d</M> is the degree of the local 
#!     minimal polynomial of <A>v</A>); 
#!  * 3rd component = the coefficients <M>a_0</M>, <M>a_1</M>, ..., 
#!     <M>a_d</M> of the local minimal polynomial; and 
#!  * 4th component = the positions of the pivots of the first component.
#! 
#! @BeginExampleSession
#! gap> A:=[ [   5,   2,  -4,   2 ],
#! >         [  -1,   0,   2,  -1 ],
#! >         [  -1,  -1,   3,  -1 ],
#! >         [ -13,  -7,  14,  -6 ] ];
#! gap> SpinMatVector(a,[1,0,0,0]);
#! [ [ [ 1, 0, 0, 0 ], [ 0, 1, -2, 1 ] ], 
#!   [ [ 1, 0, 0, 0 ], [ 5, 2, -4, 2 ] ], 
#!   [ -1, 0, 1 ],             
#!   [ 1, 2 ] ]                
#! gap> SpinMatVector(a,[0,1,0,0]);
#! [ [ [ 0, 1, 0, 0 ], [ 1, 0, -2, 1 ], [ 0, 0, 1, -1/2 ] ], 
#!   [ [ 0, 1, 0, 0 ], [ -1, 0, 2, -1 ], [ 6, 3, -4, 2 ] ], 
#!   [ 1, -1, -1, 1 ],         
#!   [ 2, 1, 3 ] ]             
#! gap> SpinMatVector(a,[1,1,0,0]);
#! [ [ [ 1, 1, 0, 0 ], [ 0, 1, 1, -1/2 ] ], 
#!   [ [ 1, 1, 0, 0 ], [ 4, 2, -2, 1 ] ], 
#!   [ 1, -2, 1 ],             
#!   [ 1, 2 ] ]                
#! @EndExampleSession
DeclareGlobalFunction("SpinMatVector");

#! @Arguments A
#! @Description
#!  Repeatedly computes the smallest invariant subspaces containing different vectors
#!  to compute a chain of cyclic subspaces. The output is a triple
#!  <C>[B,C,svec]</C>  where  <M>C</M> is such that  <M>C</M><A>A</A><M>C^-1</M>  has a block
#!  triangular shape with companion matrices along the diagonal), <M>B</M> is the
#!  row echelon form of C and svec is the list of indices where the blocks
#!  begin.
#!
#! @BeginExampleSession
#! gap> A:=[ [ 0, 1, 0, 1 ],
#! gap>      [ 0, 0, 1, 0 ],
#! gap>      [ 0, 1, 0, 1 ],
#! gap>      [ 1, 1, 1, 1 ] ];;
#! gap> sp:=CyclicChainMat(A);
#! [ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 1 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ],
#!   [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 1 ], [ 1, 1, 2, 1 ], [ 0, 0, 0, 1 ] ],
#!   [ 1, 4, 5 ] ]
#! gap> PrintArray(sp[2]*A*sp[2]^-1);
#! [ [    0,    1,    0,    0 ],  
#! [    0,    0,    1,    0 ],    
#! [    0,    3,    1,    0 ],    
#! [  1/2,  1/2,  1/2,    0 ] ]
#! @EndExampleSession
DeclareGlobalFunction("CyclicChainMat");

DeclareGlobalFunction("nfmRelMinPols");
DeclareGlobalFunction("nfmOrderPolM");
DeclareGlobalFunction("MinPolyMat");

#! @Arguments A
#! @Description
#!  Returns the minimal polynomial and a maximal vector 
#!  of the matrix <A>A</A>, that is, a vector whose local minimal polynomial 
#!  is that of <A>A</A>. This is done by repeatedly spinning up vectors until 
#!  a maximal one is found. The exact algorithm is a combination of 
#!  * the minimal polynomial algorithm by Neunhoeffer-Praeger; see <Cite Key ="Neu08"/>; and  
#!  * the method described by Bongartz 
#!      (see <Cite Key ="Bon14"/>) for computing 
#!      maximal vectors.
#!
#!  See also the article by Geck at <Cite Key ="Gec20"/>. 
#!
#! @BeginExampleSession
#! gap> A:=[ [  2,  2,  0,  1,  0,  2,  1 ],
#! >         [  0,  4,  0,  0,  0,  1,  0 ],
#! >         [  0,  1,  1,  0,  0,  1,  1 ],
#! >         [  0, -1,  0,  1,  0, -1,  0 ],
#! >         [  0, -7,  0,  0,  1, -5,  0 ],
#! >         [  0, -2,  0,  0,  0,  1,  0 ],
#! >         [  0, -1,  0,  0,  0, -1,  1 ] ];
#! gap> MaximalVectorMat(A);
#! [ [ 1, -2, 1, 1, 0, 0, 1 ], x_1^4-7*x_1^3+17*x_1^2-17*x_1+6 ]
#! gap> v:=last[1];                    
#! gap> SpinMatVector(A,v)[3];         
#! [ 6, -17, 17, -7, 1 ]                    
#! @EndExampleSession
#! In the following example, <M>M_2</M> is the (challenging) test matrix 
#! from the paper by Neunhoeffer-Praeger:
#!
#! @BeginExampleSession
#! gap> LoadPackage("AtlasRep");; g:=AtlasGroup("B",1); M2:=g.1+g.2+g.1*g.2;
#! <matrix group of size 4154781481226426191177580544000000 with 2 generators>
#! <an immutable 4370x4370 matrix over GF2>
#! gap> SetInfoLevel(Infonofoma,1);
#! gap> MaximalVectorMat(M2);;time;
#! #I Degree of minimal polynomial is 2097
#! 6725
#! gap> MinimalPolynomial(M2);;time;       
#! 13415
#! gap> LoadPackage("cvec");               
#! gap> MinimalPolynomial(CMat(M2));;time;   
#! 9721
#! @EndExampleSession
DeclareGlobalFunction("MaximalVectorMat");

#! @Arguments T,d
#! @Description
#! Modifies an already given  complementary subspace 
#! to the  complementary subspace defined by  Jacob;  concretely, this is 
#! realized by assuming that  <A>T</A>  is a matrix in block triangular shape, 
#! where the upper left diagonal block is a companion matrix (as returned 
#! by 'RatFormStep1'; the variable <A>d</A> gives the size of that block.  
#! (If <A>T</A> gives a  maximal cyclic subspace,  then  Jacob's complement is  
#! also  <A>T</A>-invariant;  but even if not,  it appears  to be  very useful 
#! because it produces many zeroes.)
DeclareGlobalFunction("JacobMatComplement");

DeclareGlobalFunction("BuildBlockDiagonalMat");
DeclareGlobalFunction("BuildBlockDiagonalMat1");

#! @Arguments A,v
#! @Description
#! Spins up a vector  <A>v</A> under a  matrix  <A>A</A>,  computes
#! a complementary subspace  (using  Jacob's construction),  and performs
#! the base change. The output is a quadruple  <C>[A1,P,pol,str]</C> where <M>A1</M> is
#! the new matrix, <M>P</M> is the base change,  <M>pol</M> is  the minimal polynomial
#! and <M>str</M> is either 'split' or 'not', according to whether the extension
#! is split or not. The second form repeatedly applies 'RatFormStep1J' in
#! order to obtain an invariant complement.
#!
#! @BeginExampleSession
#! gap> v:=[ 1, 1, 1, 1 ];;
#! gap> A:=[ [ 0, 1, 0, 1 ],
#! gap>      [ 0, 0, 1, 0 ],
#! gap>      [ 0, 1, 0, 1 ],
#! gap>      [ 1, 1, 1, 1 ] ];;
#! gap> PrintArray(RatFormStep1J(A,v)[1])
#! [ [  0,  1,  0,  0 ],    
#!   [  0,  0,  1,  0 ],    
#!   [  0,  3,  1,  0 ],    
#!   [  1,  0,  0,  0 ] ]   
#! gap> PrintArray(RatFormStep1Js(A,v)[1])";
#! [ [  0,  1,  0,  0 ],    
#!   [  0,  0,  1,  0 ],    
#!   [  0,  0,  0,  1 ],     
#!   [  0,  0,  3,  1 ] ]   
#! @EndExampleSession
DeclareGlobalFunction("RatFormStep1");
DeclareGlobalFunction("RatFormStep1J");

DeclareGlobalFunction("nfmCompanionMat");
DeclareGlobalFunction("nfmCompanionMat1");

#! @Arguments A
#! @Description
#!  Returns the invariant factors of a matrix <A>A</A>
#!  and an invertible matrix <M>P</M> such that <M>P.A.P^{{-1}}</M> is the 
#!  Frobenius normal form of <A>A</A>. The algorithm first computes a maximal 
#!  vector and an <A>A</A>-invariant complement following Jacob's construction
#!  (as described in matrix language in <Cite Key ="Gec20"/>); then the 
#!  algorithm continues recursively. It works for matrices over any field 
#!  that is available in GAP. The output is a triple with 
#!  * 1st component  = list of invariant factors; 
#!  * 2nd component = base change matrix <M>P</M>; and 
#!  * 3rd component = indices where the various blocks in the normal form 
#!       begin.
#!  You can also use  'CreateNormalForm( f[1] );' to produce the Frobenius
#!  normal form. (This function just builds the block diagonal matrix with 
#!  diagonal blocks given by the companion matrices corresponding to the 
#!  various invariant factors of <A>A</A>.) 
#! 
#! @BeginExampleSession
#! gap> A:=[ [  2,  2,  0,  1,  0,  2,  1 ],
#! >         [  0,  4,  0,  0,  0,  1,  0 ],
#! >         [  0,  1,  1,  0,  0,  1,  1 ],
#! >         [  0, -1,  0,  1,  0, -1,  0 ],
#! >         [  0, -7,  0,  0,  1, -5,  0 ],
#! >         [  0, -2,  0,  0,  0,  1,  0 ],
#! >         [  0, -1,  0,  0,  0, -1,  1 ] ];
#! gap> f:=FrobeniusNormalForm(A);
#! [ [ x_1^4-7*x_1^3+17*x_1^2-17*x_1+6, x_1^2-3*x_1+2, x_1-1 ], 
#!                                 
#!   [ [    1,   -2,    1,    1,    0,    0,    1 ],
#!     [    2,   -7,    1,    2,    0,   -1,    3 ],
#!     [    4,  -26,    1,    4,    0,   -8,    6 ],
#!     [    8,  -89,    1,    8,    0,  -35,   11 ],
#!     [ -1/2,   -2,    0,  1/2,    0,   -2, -3/2 ],
#!     [   -1,   -4,    0,    0,    0,   -4,   -2 ],
#!     [    0,  9/4,    0,   -3,    1,  5/4,  1/4 ] ],
#!                                  
#!   [ 1, 5, 7 ]  ]                 
#! gap> PrintArray(f[2]*A*f[2]^-1);
#! [ [   0,   1,   0,   0,   0,   0,   0 ], 
#!   [   0,   0,   1,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   1,   0,   0,   0 ],
#!   [  -6,  17, -17,   7,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   1,   0 ],
#!   [   0,   0,   0,   0,  -2,   3,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   1 ] ]
#! @EndExampleSession
DeclareGlobalFunction("FrobeniusNormalForm");

DeclareGlobalFunction("CreateNormalForm");
DeclareGlobalFunction("FrobeniusNormalForm1");

#! @Arguments A
#! @Description
#!  Returns the invariant factors of the matrix <A>A</A>,
#!  i.e.,  the minimal polynomials of the  diagonal blocks in the rational 
#!  canonical form  of <A>A</A>. Thus, 'InvariantFactorsMat' also specifies the
#!  rational canonical form of <A>A</A>, but without computing the base change.
#!
#! @BeginExampleSession
#! gap> InvariantFactorsMat([ [ 2,  2, 0, 1, 0,  2, 1 ],
#!                            [ 0,  4, 0, 0, 0,  1, 0 ],
#!                            [ 0,  1, 1, 0, 0,  1, 1 ],
#!                            [ 0, -1, 0, 1, 0, -1, 0 ],
#!                            [ 0, -7, 0, 0, 1, -5, 0 ],
#!                            [ 0, -2, 0, 0, 0,  1, 0 ],
#!                            [ 0, -1, 0, 0, 0, -1, 1 ] ]);
#!   #I Degree of minimal polynomial is 4
#!   #I Degree of minimal polynomial is 2
#!   [ x^4-7*x^3+17*x^2-17*x+6, x^2-3*x+2, x-1 ]
#! @EndExampleSession
DeclareGlobalFunction("InvariantFactorsMat");

DeclareGlobalFunction("nfmFrobInv");
DeclareGlobalFunction("SquareFreePol");

#! @Arguments A,f
#! @Description
#!  Returns the unique pair of matrices <M>D</M>, 
#!  <M>N</M> such that the matrix <A>A</A> is written as <M>A=D+N</M>, where 
#!  <M>N</M> is a nilpotent matrix and <M>D</M> is a matrix that is 
#!  diagonalisable (over some extension field of the default field of 
#!  <A>A</A>), such that <M>D.N=N.D</M>; the argument <A>f</A> is a 
#!  polynomial such that <M>f(A)=0</M> (e.g., the minimal polynomial of 
#!  <A>A</A>). This is called the Jordan-Chevalley decomposition of <A>A</A>; 
#!  the algorithm is based on <Cite Key ="Gec22"/>. Note that this 
#!  algorithm does not require the knowledge of the eigenvalues of <A>A</A>; 
#!  it works over any perfect field that is available in GAP.
#!
#! @BeginExampleSession
#! gap> A:=[ [  6, -2,  6,  1,  1 ],
#! >         [  1, -1,  2,  1, -2 ],
#! >         [ -2,  0, -1,  0, -1 ],
#! >         [ -1,  0, -2,  2, -1 ],
#! >         [ -4,  4, -6, -2,  3 ] ];
#! gap> jc:=JordanChevalleyDecMat(A,MinimalPolynomial(A));
#! [ [ [  4,  0,  4, -1,  1 ], 
#!     [  1,  0,  1,  1, -1 ], 
#!     [ -1, -1,  0,  1, -1 ], 
#!     [  0,  0, -2,  3,  0 ], 
#!     [ -3,  2, -4, -1,  2 ] ], 
#!   [ [  2, -2,  2,  2,  0 ], 
#!     [  0, -1,  1,  0, -1 ], 
#!     [ -1,  1, -1, -1,  0 ], 
#!     [ -1,  0,  0, -1, -1 ], 
#!     [ -1,  2, -2, -1,  1 ] ] ]
#! gap> MinimalPolynomial(jc[1]);
#! x_1^3-5*x_1^2+9*x_1-5
#! gap> Factors(last);
#! [ x_1-1, x_1^2-4*x_1+5 ]  
#! gap> MinimalPolynomial(jc[2]);
#! x_1^2                     
#! @EndExampleSession
#!  If the input matrix is very large, then 'JordanChevalleyDecMatF(<A>A</A>);' 
#!  may be more efficient; this function first computes the Frobenius normal 
#!  form of <A>A</A> and then applies 'JordanChevalleyDecMat' to each diagonal 
#!  block. (The result will be the same as that of 
#!  'JordanChevalleyDecMat(<A>A</A>);)'
DeclareGlobalFunction("JordanChevalleyDecMat");

#! @Arguments A
#! @Description
#!  First computes the Frobenius normal form and
#!  then applies 'JordanChevalleyDecMat' to each diagonal block.
DeclareGlobalFunction("JordanChevalleyDecMatF");

DeclareGlobalFunction("CheckFrobForm");
DeclareGlobalFunction("CheckJordanChev");

DeclareGlobalFunction("nfmmat1");

#! @Section Further documentation
#! The above functions, as well as a number of further auxiliary functions, 
#! are all contained and defined in the file 'pgk/nofoma-1.0/gap/nofoma.gi'; 
#! in that file, you can also find further inline documentation for the 
#! auxiliary functions.



