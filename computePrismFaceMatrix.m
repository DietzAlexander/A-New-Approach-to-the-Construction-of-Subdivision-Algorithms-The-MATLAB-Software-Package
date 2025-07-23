function [FaceMatrix] = computePrismFaceMatrix(n)
%COMPUTEPRISMFACEMATRIX Computes a face matrix for an n-sided prism. The
%sides of the prism have valence 4 where the top and the bottom of the
%prism have valence n
%
%Input: n: Amount of sides of the prism (without top and bottom)
%Output: Face matrix for the n-sided prism. 
%
%computeFaceMatrix(n)  Computes a face matrix for an n-sided prism.

%Input check. Depending on the input the return variable needs another size
if n<2
    error("A prism cannot have " + n + " side(s).")
elseif n==2
    warning('A prism with 2 sides might not be well defined.')

    %Biggest face has 4 vertices (side faces)
    FaceMatrix=zeros(n+2,4);
elseif n==3
    %Biggest face has 4 vertices (side faces)
    FaceMatrix=zeros(n+2,4);
else
    %Biggest face has n vertices (top and bottom face)
    FaceMatrix=zeros(n+2,n);
end

%Sets the top face. Numbers it from 1 to n.
FaceMatrix(1,1:n)=1:n;

%Sets the bottom face. Numbers it from n+1 to 2n.
FaceMatrix(2,1:n)=n+1:2*n;

%Sets the side faces
for i=1:n
    %Every side face has 2 points of the top face and 2 points of the
    %bottom face. The numeration is [i, i+1, n+i+1, n+i]
   FaceMatrix(i+2,1:4)=[i,mod(i,n)+1,n+mod(i,n)+1,n+i]; 
end

end


