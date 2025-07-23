function [ S ] = computeBiQuadraticSubdivisionMatrixV1Big( n )
%computeBiQuadraticSubdivisionMatrixV1Big computes the subdivision matrix
%for the quadratic 2-D case with the size of a  ring.
%
% Input:  n: Amount of points of the central 2-polytope
% Output: S: 9n x 9n subdivision matrix
%
%computeBiQuadraticSubdivisionMatrixV1Big(n) computes the subdivision matrix for 
% the cubic 2-D case with the size of a ring.

    %The Amount of control points (9 per point in the center)
    AmountOfPoints = 9*n;
    
    %the empty subdivision matrix
    S=zeros(AmountOfPoints);
  
    %the center part of the big subdivision matrix
    SCenter=computeBiQuadraticSubdivisionMatrixV1(n);
    
    %The main ounter block
    S(1:8,1:9)=...
    [0,          0,          0,          0,          9/16,       3/16,       0,              3/16,           1/16;...
     0,          0,          0,          0,          3/16,       9/16,       0,              1/16,           3/16;...
     0,          0,          0,          0,          0,          9/16,       0,              0,              3/16;...
     0,          0,          0,          0,          3/16,       1/16,       0,              9/16,           3/16;...
     0,          0,          0,          0,          1/16,       3/16,       0,              3/16,           9/16;...
     0,          0,          0,          0,          0,          3/16,       0,              0,              9/16;...
     0,          0,          0,          0,          0,          0,          0,              9/16,           3/16;...
     0,          0,          0,          0,          0,          0,          0,              3/16,           9/16];

    %side part of the outer block
    S(1:8,17:18)=...
    [0,          0;...
     0,          0;...
     3/16,       1/16;...
     0,          0;...
     0,          0;...
     1/16,       3/16;...
     0,          0;...
     0,          0];

     %additional side part of the outer block
    S(1:8,end-3)=[0;  0;  0;  0;  0;    0;  3/16;   1/16];
    S(1:8,end)=  [0;  0;  0;  0;  0;    0;  1/16;   3/16];
    
    for j=1:n   
           S(9,(j-1)*9+9)=SCenter(1,j); 
    end

    
    
    
    
    %Due to symmetry the created part just has to be shifted
    for i=2:n
        S(9*(i-1)+1:9*(i),:)=circshift(S(9*(i-2)+1:9*(i-1),:),9,2);
    end


end

