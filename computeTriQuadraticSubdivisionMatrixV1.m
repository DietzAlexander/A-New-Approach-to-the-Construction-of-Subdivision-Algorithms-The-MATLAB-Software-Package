function [S] = computeTriQuadraticSubdivisionMatrixV1(Input,mu,varargin)
%COMPUTETRUQUADRATICSUBDIVISIONMATRIXV1 computes out of a polytope structure the
%corresponding subdivision matrix, which refines a polytope to a smaller one.
%
% Input:    adjacency matrix or face matrix or edge matrix of the
%           polytope
%
% mu:       The subsubdominant eigenvalue
%
% varargin: Several optional commands (see the manual)
%
% S:        Subdivision matrix
%
% For more information see the corresponding manual

%Number of additional inputs
NumberOfAdditionalInput = nargin - 2;

%Standard values for several options of varargin
status=false;

for i=1:NumberOfAdditionalInput
    if ischar(varargin{i}) && strcmp(varargin{i},'Status')
        status=true;
    end
end

%Check, if mu is numeric
if ~isnumeric(mu)
    error('mu is not a numeric value.')
end


%Checks whether the input is a face matrix, an adjacency matrix or an edge matrix 
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


%compute the 3-polytope
[PointCoordinates3D,~,StatusString,statusCounter] = Construct3Polytope(Input,varargin{:});


if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Check, if the input is a prism.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%-------------------------------------------------------------------------
% Checks, if the input is a prism or not
%-------------------------------------------------------------------------

%if the amount of vertices is odd, it is not a prism
if mod(AmountOfVertices,2)==1 || AmountOfVertices<6
    IsPrismWahrheitswert=false;
else
    %creates a graph for the input 
    AA=graph(AdjacencyMatrix);

    %and for a prism of same size
    BB=graph(FaceToAdjacencyMatrix(computePrismFaceMatrix(AmountOfVertices/2)));

    %if there is an isomorphism, it is a a prism
    I=isomorphism(AA,BB);
    if ~isempty(I)
        IsPrismWahrheitswert=true;
    else
        IsPrismWahrheitswert=false;
    end
end

if IsPrismWahrheitswert
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Input is a prism.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end
else
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Input is not prism.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end
end

%-------------------------------------------------------------------------
% Computes the subdivision matrix
%-------------------------------------------------------------------------

%if it is a prism...
if IsPrismWahrheitswert

    %compute the standard Doo-Sabin subdivision matrix
    S=computeBiQuadraticSubdivisionMatrixV1(AmountOfVertices/2);

    %compute the Kroneckerproduct
    S=1/4*[3*S,S;S,3*S];
    S=S(I,I);
else

    %if it is not a prism get the primal points
    PointCoordinates3D=PointCoordinates3D(1:AmountOfVertices,:);

    %center them, that the barycenter of th primal points is 0
    PointCoordinates3D=PointCoordinates3D-sum(PointCoordinates3D)/length(PointCoordinates3D);
    
    %get a orthogonal version of the points
    PointCoordinates3D=orth(PointCoordinates3D);
    
    %compute the subdivision matrix
    S=(ones(AmountOfVertices,1)*ones(AmountOfVertices,1)'/norm(ones(AmountOfVertices,1))^2+(1-2*mu)/(2-2*mu)*PointCoordinates3D*(PointCoordinates3D'))*(1-mu)+eye(AmountOfVertices)*mu;
end

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Subdivision matrix was computed.\n'];
    fprintf(StatusString);
end


end

