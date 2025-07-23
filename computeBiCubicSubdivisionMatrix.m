function [S] = computeBiCubicSubdivisionMatrix(n)
%computeBiCubicSubdivisionMatrix computes the subdivision matrix for the cubic
% 2-D case.
%
% Input:  n: Amount of points of the dual 2-polytope
% Output: S: 2n+1 x 2n+1 subdivision matrix
%
%computeDooSabinSubdivisionMatrix(n) computes the subdivision matrix for 
% the cubic 2-D case.

%----------------------------------------------------------------------
% Compute the eigenstructure
%----------------------------------------------------------------------

%Compute the tangent points
r=zeros(1,n);
PTangent=[real(exp(2*pi*1i*((0:n-1)+r)/n))',imag(exp(2*pi*1i*((0:n-1)+r)/n))'];

GLS=zeros(n,2);
RS=zeros(n,1);

%Prepare the GLS for the corner points
for i=1:n
    GLS(i,:)=PTangent(i,:);
    RS(i)=PTangent(i,:)*PTangent(i,:)';
end

CornerPoints=zeros(n,2);

%Solve the GLS to compute the corner points
for i=1:n
    if i<n
        CornerPoints(i,:)=GLS([i,i+1],:)\RS([i,i+1]);
    else
        CornerPoints(i,:)=GLS([n,1],:)\RS([n,1]);
    end
    
end

%Add the central point.
PointCoordinates=[CornerPoints;PTangent;0,0];

%----------------------------------------------------------------------
% Compute the matrix MBar
%----------------------------------------------------------------------
MBar=zeros(4*n,2*n+1);
II=zeros(4*n,2);
for i=1:n

    %Gets the points of each square and its indices
    if i==n
        PLocal=PointCoordinates([i,2*n,2*n+1,1+n],:)-1/2*PointCoordinates(i,:);
        I=[i,2*n,2*n+1,1+n]';
    else
        PLocal=PointCoordinates([i,i+n,2*n+1,i+n+1],:)-1/2*PointCoordinates(i,:);
        I=[i,i+n,2*n+1,i+n+1]';
    end

    %Prepare input for the algorithm which computes the CDVMatrix
    Polygon=create2PolytopeAdjacencyMatrix(4);
    Points=PLocal;
    Center=[0,0];
    
    %Compute CDV Matrix for the cube
    [M]=computeCDVMatrix2D(Polygon,Points,Center);

    %Set the values
    M=diag(1./sum(M,2))*M;
    MBar(4*(i-1)+1:4*i,I)=M;
    II(4*(i-1)+1:4*i,:)=[I,[i;i;i;i]];
end

%----------------------------------------------------------------------
% Compute the matrix BBar
%----------------------------------------------------------------------

BBar=zeros(2*n+1,4*n);
for i=1:2*n+1
    I=find(II(:,1)==i);

    if i<=n
        BBar(i,I)=1;
    elseif i<=2*n
        BBar(i,I)=1/2;
    else
        BBar(i,I)=1/n;
    end
end


%----------------------------------------------------------------------
% Compute the matrix C and S
%----------------------------------------------------------------------

C=BBar*MBar;

C=2*C;
L=zeros(2*n+1);
R=zeros(2*n+1);
for i=1:2*n+1
    for j=1:2*n+1
        if i>j
            L(i,j)=C(i,j);
        elseif i<j
            R(i,j)=C(i,j);
        end
    end
end

for i=1:2*n+1
    L(i,i)=1-sum(L(i,:));
    R(i,i)=1-sum(R(i,:));
end

L=expm(log(2)*(L-eye(2*n+1)));
R=expm(log(2)*(R-eye(2*n+1)));

S=L*R;


end






%Also a local function in computeTriCubicSubdivisionMatrix
function [CDVNormalized,CDV] = computeCDVMatrix2D(AdjacencyMatrix,Controlpoints,Center)
%COMPUTECDVMATRIX2D  computes a 2-dimensional CDV matrix for the given
%controlpints, whichs 2-polytope lies on a circle with mitpoint Center.
%
%Input:     AdjacencyMatrix:    The adjacency matrix of the structure
%           Controlpoints:      The controlpoints of the structure
%           Center:             The midpoint of the circle on which the
%                               2-polytope lies
% Output:   CDVNormalized:      The normalized CDV matrix
%           CDV:                The CDV matrix


%The amount of the adjacency matrix
[a,~]=size(AdjacencyMatrix);

%Initializing the CDV matrix
CDV=zeros(a);

%the size of the control points
[c,d]=size(Controlpoints);

%as this algorithm wors in R^3...
if d==2
    Controlpoints=[Controlpoints,zeros(c,1)];
end

%as this algorithm wors in R^3...
[c,d]=size(Center);
if d==2
    Center=[Center,zeros(c,1)];
end

%Shifts the control points such that the center is in [0,0,0]
PShifted=Controlpoints-Center;

%computes each entry
for i=1:a
    for j=1:a

        %entries are just unequal to 0, if the adjacency matrix is unequal
        %to 0
        if AdjacencyMatrix(i,j)==1

           %compute the entries
           Ui=[PShifted(i,:)];
           Uj=[PShifted(j,:)];
           UU=cross(Ui,Uj);
           CDV(i,j)=-1/(sqrt(UU*UU'));
  
        end
    end
    %unnecessary condition to avoid dividing by 0
    if sum(AdjacencyMatrix(i,:))>0

        %compute the duaginal entries
        UStrich=CDV(i,:)*PShifted;
        CDV(i,i)=-UStrich/PShifted(i,:);
    end

end

%normalize the CDV matrix
CDVNormalized=CDV./sum(CDV,2);


end




%Also a local function in computeTriCubicSubdivisionMatrix
function [AdjacencyMatrix] = create2PolytopeAdjacencyMatrix(n)
% CREATE2POLYTOPEADJACENCYMATRIX creates an adjacency matrix for a
% 2-polytope with n vertices
%
% Input:    n:                  amount of vertices
% Output:   AdjacencyMatrix:    the adjacency matrix
%
% create2PolytopeAdjacencyMatrix(n) creates an adjacency matrix for a
% 2-polytope with n vertices

%initializing the adjacency matrix
AdjacencyMatrix=zeros(n);

%goes through each entry
for i=1:n
    for j=1:n

        %if one of the two conditions is fulfilled the entry is 1
        if j==i+1
            AdjacencyMatrix(i,j)=1;
        end
        if j==i-1
            AdjacencyMatrix(i,j)=1;
        end
    end
end

%the two special cases (mod n)
AdjacencyMatrix(1,n)=1;
AdjacencyMatrix(n,1)=1;
end