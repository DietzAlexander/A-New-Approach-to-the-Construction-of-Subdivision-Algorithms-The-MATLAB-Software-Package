function [DualAdjacency,PolygonFound] = CreateTestScenario(n,maxVertices)
%CREATETESTSCENARIO computes an adjacency matrix of a randomly generated
%planar 3-connected graph.
%
%Input: n:          The maximum amount of faces of the graph
%       maxVertices The maximum amount of vertices of the graph
%
%Output: DualAdjacency: The  adjacency matrix of the randomly genarated  
%                       graph if the algorithm was sucessfull and [] if
%                       not.
%        PolygonFound:  is true, if the algorithm was sucessfull and false,
%                       if not.
%
%CreateTestScenario(n,maxVertices) computes an adjacency matrix of a 
%randomly generated planar 3-connected graph.



%-----------------------------------------------------------------
%   Setting start values
%-----------------------------------------------------------------

%Random generator is based on current time
rng("shuffle");

PolygonFound=false;
counter=1;

DualAdjacency=[];


%Choose randomly, if edges of the dual polytope are deleted. If there is no
%deletion, the final polytope belongs to a hexahedral structure (all nodes 
%have exact 3 neighbors).
r=randi(2);
if r==2
    EdgeDeletionNumber=n;
else
    EdgeDeletionNumber=0;
end

%The whole loop runs as long as a random polytope was found or as long as 
% 20 try were done 
while ~PolygonFound && counter<20

    %-----------------------------------------------------------------
    %   Compute Random Points on the shpere
    %-----------------------------------------------------------------
    
    %Variable that determines, if a set of points on the sphere was found or not
    PointsOnSphereFound=false;
    
    while ~PointsOnSphereFound
    
        %Create a random set of Points. Each value of each coordinate is
        %between -1 and 1
        ControlPointsTriangulation=2*(rand(n,3)-0.5);
    
        %Try to normalize the points. This is possible, if there is no point accidentally on the origin. 
        if min(abs(sum(ControlPointsTriangulation,2)))>0
    
            %If normalisation is possible normalize
            ControlPointsTriangulation=ControlPointsTriangulation./ vecnorm(ControlPointsTriangulation,2,2);
    
            %and set the found parameter to true
            PointsOnSphereFound=true;
        end
    end
    
    %-----------------------------------------------------------------
    %   Compute a delaunay triangulation of these points
    %-----------------------------------------------------------------
    
    
    %Compute the delaunay triangulation
    DT=delaunayTriangulation(ControlPointsTriangulation);
        
    %Compute the convex hull of the delaunay triangulation, which is a face
    % matrix
    [K,~] = convexHull(DT);
    
    %compute the adjacency matrix out of the face matrix
    AdjacencyMatrix=FaceToAdjacencyMatrix(K);

    %delete all empty columns of the adjacency matrix
    AdjacencyMatrix(sum(AdjacencyMatrix,2)==0,:)=[];

    %delete all empty rows of the adjacency matrix
    AdjacencyMatrix(:,sum(AdjacencyMatrix)==0)=[];

    %get the amount of vertices of the graph
    AmountOfVertices=length(AdjacencyMatrix);

    %get the amount of vertices of the graph
    AmountOfEdges=sum(sum(AdjacencyMatrix))/2;
    
    %decide, how many edges of the graph at most should be deleted (the 
    %minimal graph, which is a tetraedron has 6 edges)
    CurrentDeletionNumber=min(AmountOfEdges-6,EdgeDeletionNumber);
    
    %the amount of maximum trials of edge deletion
    MaxTrials=CurrentDeletionNumber^2;

    %tries to delete an edge as long as the maximum of edges are deleted or
    %as log as the MaxTrials are udes up
    while CurrentDeletionNumber >0 && MaxTrials>0
        i=randi(AmountOfVertices);
        j=randi(AmountOfVertices);
        if AdjacencyMatrix(i,j)==1
            AdjacencyMatrix(i,j)=0;
            AdjacencyMatrix(j,i)=0;
            CurrentDeletionNumber=CurrentDeletionNumber-1;
        end
        MaxTrials=MaxTrials-1;
    end
    
    %-------------------------------------------------------------------------
    %       Input Check 
    %-------------------------------------------------------------------------
    
    %checks several properties of the generated AdjacencyMatrix
    if length(K)>maxVertices
    
    elseif sum(diag(AdjacencyMatrix)==zeros(AmountOfVertices,1))~=AmountOfVertices
 
    elseif AmountOfVertices<4

    elseif sum(sum(AdjacencyMatrix))<8

    elseif ~is3Connected(AdjacencyMatrix)

    elseif ~PlanarityTest(AdjacencyMatrix)

    
    %If all the properties are fulfilled, once can generate the dual
    %polytope
    else

        %-------------------------------------------------------------------------
        %       Constructing the dual graph
        %-------------------------------------------------------------------------

        PrimalGraph=graph(AdjacencyMatrix);
        
        %Polyeder theorem of Euler
        PrimalEdges=PrimalGraph.Edges.EndNodes;
        
        %List of Faces of the Polytope
        FaceMatrix=computeFaceMatrix(AdjacencyMatrix,false);
        

        %Building the dual graph and the quad graph
        DualGraph=graph([]);
        for i=1:length(PrimalEdges)

            %Find the 2 faces, that have one edge in common.
            [row1,~]=find(FaceMatrix==PrimalEdges(i,1));
            [row2,~]=find(FaceMatrix==PrimalEdges(i,2));
            Connection=row1(ismember(row1,row2));

            %Add a line between these 2 faces which are the nodes in the dual graph
            DualGraph=addedge(DualGraph,Connection(1),Connection(2),0);
        end

        DualAdjacency=DualGraph.adjacency; 

        DualAmountOfVertices=length(DualAdjacency);

        %Check, if the dual graph fulfill all properties
        if sum(diag(DualAdjacency)==zeros(DualAmountOfVertices,1))~=DualAmountOfVertices
        
        elseif DualAmountOfVertices<4
            
        elseif sum(sum(DualAdjacency))<8
            
        elseif ~is3Connected(DualAdjacency)
            
        elseif ~PlanarityTest(DualAdjacency)
            
        else
            %if yes, the polygon was found.
            PolygonFound=true;
        end
        

    end
    EdgeDeletionNumber=EdgeDeletionNumber-1;
    counter=counter+1;

end

%if no polytope could be found, the variable DualAdjacency is erased
if ~PolygonFound
    DualAdjacency=[];
end
