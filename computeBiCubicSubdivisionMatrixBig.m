function [ S ] = computeBiCubicSubdivisionMatrixBig(n)
%computeBiCubicSubdivisionMatrix computes the subdivision matrix for the cubic
% 2-D case with the size of a double ring.
%
% Input:  n: Amount of points of the central dual 2-polytope
% Output: S: 30n+1 x 30n+1 subdivision matrix
%
%computeBiCubicSubdivisionMatrixBig(n) computes the subdivision matrix for 
% the cubic 2-D case with the size of a double ring.


%THe Amount of Control Points
AmountOfPoints = 30*n + 1;

%The subdivision matrix
S= zeros(AmountOfPoints);

%regular stencils
M1=[1/4,1/4,1/4,1/4];
M2=[1/16,3/8,1/16,1/16,3/8,1/16];
M22=[1/16,1/16,3/8,3/8,1/16,1/16];
M3=[1/64,3/32,1/64,3/32,9/16,3/32,1/64,3/32,1/64];

%goes through each segment
for i=1:n

    %intern number where to start in the matrix
    s=30*(i-1);

    %Fills the regular part of the matrix with the stencils
    S(2+s,[16,17,22,23]+s)=M1;
    S(3+s,[16,17,18,22,23,24]+s)=M2;
    S(4+s,[17,18,23,24]+s)=M1;
    S(5+s,[17,18,19,23,24,25]+s)=M2;
    S(6+s,[18,19,24,25]+s)=M1;
    if i~=n
        S(7+s,[18,19,58,24,25,59]+s)=M2;
    else
        S(7+s,[18,19,28-s,24,25,29-s]+s)=M2;
    end

    S(8+s,[16,17,22,23,28,29]+s)=M22;
    S(9+s,[16,17,18,22,23,24,28,29,30]+s)=M3;
    S(10+s,[17,18,23,24,29,30]+s)=M22;
    S(11+s,[17,18,19,23,24,25,29,30,31]+s)=M3;
    S(12+s,[18,19,24,25,30,31]+s)=M22;
    if i~=n
        S(13+s,[18,19,58,24,25,59,30,31,60]+s)=M3;
    else
        S(13+s,[18,19,28-s,24,25,29-s,30,31,30-s]+s)=M3;
    end

    S(14+s,[22,23,28,29]+s)=M1;
    S(15+s,[22,23,24,28,29,30]+s)=M2;
    S(16+s,[23,24,29,30]+s)=M1;
    S(17+s,[23,24,25,29,30,31]+s)=M2;
    S(18+s,[24,25,30,31]+s)=M1;
    if i~=n
        S(19+s,[24,25,59,30,31,60]+s)=M2;
    else
        S(19+s,[24,25,29-s,30,31,30-s]+s)=M2;
    end

    if i~=1
        S(20+s,[22,23,28,29,19-30,25-30]+s)=M22;
        S(21+s,[22,23,24,28,29,30,19-30,25-30,31-30]+s)=M3;
        S(22+s,[23,24,29,30,25-30,31-30]+s)=M22;
        S(23+s,[23,24,25,29,30,31,25-30,31-30,1-s]+s)=M3;
        S(24+s,[24,25,30,31,31-30,1-s]+s)=M22;
        if i~=n
            S(25+s,[24,25,59,30,31,60,31-30,1-s,61]+s)=M3;
        else
            S(25+s,[24,25,29-s,30,31,30-s,31-30,1-s,31-s]+s)=M3;
        end
    else
        S(20+s,[22,23,28,29,end-12-s,end-6-s]+s)=M22;
        S(21+s,[22,23,24,28,29,30,end-12,end-6-s,end-s]+s)=M3;
        S(22+s,[23,24,29,30,end-6-s,end-s]+s)=M22;
        S(23+s,[23,24,25,29,30,31,end-6-s,end-s,1-s]+s)=M3;
        S(24+s,[24,25,30,31,end-s,1-s]+s)=M22;
        if i~=n
            S(25+s,[24,25,59,30,31,60,end-s,1-s,61]+s)=M3;
        else
            S(25+s,[24,25,29-s,30,31,30-s,end-s,1-s,31-s]+s)=M3;
        end
    end


    if i~=1
        S(26+s,[28,29,19-30,25-30]+s)=M1;
        S(27+s,[28,29,30,19-30,25-30,31-30]+s)=M2;
        S(28+s,[29,30,25-30,31-30]+s)=M1;
        S(29+s,[29,30,31,25-30,31-30,1-s]+s)=M2;
        
    else
        S(26+s,[28,29,end-12-s,end-6-s]+s)=M1;
        S(27+s,[28,29,30,end-12,end-6-s,end-s]+s)=M2;
        S(28+s,[29,30,end-6-s,end-s]+s)=M1;
        S(29+s,[29,30,31,end-6-s,end-s,1-s]+s)=M2;
        
    end


end

%Compute the central part of the subdivision matrix
SCentral=computeBiCubicSubdivisionMatrix(n);

%Adjust the counting of the central part.
V=zeros(2*n+1,1);

V(1)=2*n+1;
counter=2;
for i=n:-1:1
    V(counter)=i;
    counter=counter+1;
    V(counter)=n+i;
    counter=counter+1;
end

%Get the coresponding counting in the Big matrix
VBig=zeros(2*n+1,1);

 VBig(1)=1;
 for i=1:n
     VBig(1+2*(i-1)+1)=30+30*(i-1);
     VBig(1+2*(i-1)+2)=31+30*(i-1);
 end

%Adjust the central matrix 
SCentral=SCentral(V,V);

%set the values
S(VBig,VBig)=SCentral;

end


