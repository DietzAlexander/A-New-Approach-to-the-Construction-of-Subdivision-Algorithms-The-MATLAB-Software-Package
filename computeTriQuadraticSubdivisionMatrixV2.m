function [S] = computeTriQuadraticSubdivisionMatrixV2(Input,varargin)
%COMPUTETRUQUADRATICSUBDIVISIONMATRIXV2 computes out of a polytope structure the
%corresponding subdivision matrix, which refines a polytope to a smaller one.
%
% Input:    adjacency matrix or face matrix or edge matrix of the
%           polytope
%
% varargin: Several optional commands (see the manual)
%
% S:        Subdivision matrix
%
% For more information see the corresponding manual

%Number of additional inputs
NumberOfAdditionalInput = nargin - 1;

%Standard values for several options of varargin
status=false;

for i=1:NumberOfAdditionalInput
    if ischar(varargin{i}) && strcmp(varargin{i},'Status')
        status=true;
    end
end

%Checks whether the input is a face matrix, an adjacancy matrix or an edge matrix 
%also check some basic plausibilitys if the input is valid
%creates then the adjacency matrix for all input types

[rows,cols]=size(Input);

%If the input has just 2 columns then we have an edge matrix
if cols ==2

    EdgeMatrix=Input;

    %The amount of vertices
    AmountOfVertices=max(max(EdgeMatrix));
    
    %Comparing the labeling of the vertices with 1,2,...,AmountOfVertices
    if sum((1:AmountOfVertices)==(unique(EdgeMatrix))') ~= AmountOfVertices
        error('Vertices are not named uniformly.');
    end

    G=graph(EdgeMatrix(:,1),EdgeMatrix(:,2));
    AdjacencyMatrix=full(G.adjacency);


%If the Input is square and every entry is either 0 or 1, then we have an adjacency matrix 
elseif rows==cols && sum(sum((Input==0)+(Input==1)))==rows^2
    AdjacencyMatrix=Input;
    
    %The amount of vertices
    AmountOfVertices=length(AdjacencyMatrix);

    if sum(sum(AdjacencyMatrix==AdjacencyMatrix'))~=rows^2
        error('Adjacency matrix as input must be symmetric.')
    end

%Otherwise we have a face matrix (or something invalid)
else
    %The amount of vertices
    AmountOfVertices=max(max(Input));
    
    %Comparing the labeling of the vertices with 1,2,...,AmountOfVertices
    if sum((1:AmountOfVertices)==(unique(Input(Input>0)))')~=AmountOfVertices
        error('Vertices are not named uniformly.');
    end

    
    AdjacencyMatrix=FaceToAdjacencyMatrix(Input);

end


PrimalGraph=graph(AdjacencyMatrix);


%Polyeder theorem of Euler
PrimalEdges=PrimalGraph.Edges.EndNodes;
AmountOfEdges=length(PrimalEdges);

[PointCoordinates3D,KiteAll,StatusString,statusCounter] = Construct3Polytope(Input,varargin{:});

%-------------------------------------------------------------------------
%       Compute the Colin de Verdiere Matrix
%-------------------------------------------------------------------------

%The Colin de Verdiere Matrix
C=zeros(AmountOfVertices);

%Compute for every edge the entry of the CDV Matrix
for i=1:AmountOfEdges
    %Get the point numbers for the primal points
    currentEdge=PrimalEdges(i,:);

    %Get the point numbers for the edge
    currentEdgePointNumber=AmountOfVertices+i;

    %Get the point numbers for the dual points
    [KitesWithThisEdge,~]=find(KiteAll==currentEdgePointNumber);
    DualPoints=unique(KiteAll(KitesWithThisEdge,3));

    %A just teoretical error. Shount not be thrown at any point but if its
    %thrown something with the code is not ok
    if length(DualPoints)~=2
        error('Something went wrong with the kite structure. One edgepoint has more or less then 2 dual points.')
    end

    %The four points needed for calculation
    ui=PointCoordinates3D(currentEdge(1),:);
    uj=PointCoordinates3D(currentEdge(2),:);
    wf=PointCoordinates3D(DualPoints(1),:);
    wg=PointCoordinates3D(DualPoints(2),:);

    %Calculation of the matrix entry
    MatrixEntry=-abs(norm(wf-wg))/abs(norm(cross(ui,uj)));
    

    %Setting the matrix entries
    C(currentEdge(1),currentEdge(2))=MatrixEntry;
    C(currentEdge(2),currentEdge(1))=MatrixEntry;
end

%Compute the sum of the off diagonal entries multiplied with the
%coordiantes
SumMPointsM=C*PointCoordinates3D(1:AmountOfVertices,:);

%compute the coordinates of the diagonal entries
for i=1:AmountOfVertices
    C(i,i)=-SumMPointsM(i,:)/PointCoordinates3D(i,:);
end

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Colin de Verdiere matrix was computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%-------------------------------------------------------------------------
%       Adjust the eigenspectrum
%-------------------------------------------------------------------------

%Compute the row sum of the Colin de Verdiere Matrix C
N=diag(1./sum(C,2));

%Norm the row of M
TildeC=N*C;

%Compute the final subdivision matrix as exponantial of TildeM
S=expm(log(2)*(TildeC-eye(AmountOfVertices)));

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Spectrum was adjusted.\n'];
    fprintf(StatusString);
end

end