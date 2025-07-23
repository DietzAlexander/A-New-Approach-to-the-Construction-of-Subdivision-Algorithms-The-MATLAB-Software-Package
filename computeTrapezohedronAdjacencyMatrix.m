function [A] = computeTrapezohedronAdjacencyMatrix(n)
% COMPUTETRAPEZOHEDRONADJACENCYMATRIX computes the adjacency matrix for a
% trapezohedron with 2n+2 vertices.
%
% Input:    Amount of vertices: n -> 2*n+2
% Output:   Adjacency matrix of a trapezohedron
%
% computeTrapezohedronAdjacencyMatrix(n) computes the adjacency matrix for a
% trapezohedron with 2n+2 vertices.

% The central part of the adjacency matrix
T2n=zeros(2*n);

%fills the central part
for i=1:(2*n-1)
    T2n(i,i+1)=1;
end

%fill the central part
for i=2:(2*n)
    T2n(i,i-1)=1;
end

%fill the central part
T2n(1,2*n)=1;
T2n(2*n,1)=1;

%fill the upper strip
StripTop=zeros(2*n,1);
StripTop(1:2:2*n)=1;

%fill the bottom strip
StripBottom=zeros(2*n,1);
StripBottom(2:2:2*n)=1;

%compute the adjacency matrix
A=[0,StripTop',0;...
    StripTop,T2n,StripBottom;...
    0,StripBottom',0];
end
