function [Result] = is3Connected(AdjacencyMatrix)
%IS3CONNECTED Checks, if the graph related to the AdjacencyMatrix is 
%3-connected or not. Returns true if 3-connected and false, if not.
%
%Input: n x n adjacency matrix
%Output: Boolean true or false
%
%is3Connected(AdjacencyMatrix) checks if the graph to AdjacencyMatrix is
%3-connected.

%The amount of vertices
AmountOfVertices=length(AdjacencyMatrix);

%The minimal value a vertex is connected to another vertex. At the
%beginning the value is AmountOfVertices-1 the value for a complete graph.
%This value might be reduced during the algorithm to the minimal
%connectivity
MinConnectivity=AmountOfVertices-1;

%Size of the adjacency matrix to check whether it is square
[rows,cols]=size(AdjacencyMatrix);


%-------------------------------------------------------------------
%  Input Check
%-------------------------------------------------------------------

%Square check of the adjacency matrix
if rows ~= cols
    error('AdjacencyMatrix has to be square.');
end

%Check for self loops
if sum(diag(AdjacencyMatrix)==zeros(AmountOfVertices,1))~=AmountOfVertices
    error('Self loops are permited.')
end

%Check for symmetry
if sum(sum(AdjacencyMatrix==AdjacencyMatrix'))~=rows^2
    error('AdjacencyMatrix as input must be symmetric.')
end

%-------------------------------------------------------------------
%  Construction of the directed graph
%-------------------------------------------------------------------

%Variable of the extended graph
ExtendedAdjacencyMatrix=zeros(2*AmountOfVertices);


%Replaces every vertex v by 2 vertices v' v'' such that
% (vertices linking to v)->v'->v''->(vertices v links to)
for k=1:AmountOfVertices
    ExtendedAdjacencyMatrix(k,k+AmountOfVertices)=1;
    for l=1:AmountOfVertices
        if AdjacencyMatrix(k,l)==1
            ExtendedAdjacencyMatrix(k+AmountOfVertices,l)=inf;
            ExtendedAdjacencyMatrix(l+AmountOfVertices,k)=inf;
        end
    end
end

%creates a data structure of directed edges form Edges(:,1) to Edges(:,2)
% with weight Edges(:,3) which is 1 for the connection between the doubled
% vertices and infinity on each other edge
AmountOfEdges=nnz(ExtendedAdjacencyMatrix);
Edges=zeros(AmountOfEdges,3);
counter=1;
for k=1:2*AmountOfVertices    
    for l=1:2*AmountOfVertices
        if ExtendedAdjacencyMatrix(k,l)~=0
            Edges(counter,1)=k;
            Edges(counter,2)=l;
            Edges(counter,3)=ExtendedAdjacencyMatrix(k,l);
            counter=counter+1;
        end
    end
end

%Make a graph out of the Edges
G=digraph(Edges(:,1),Edges(:,2),Edges(:,3));

%-------------------------------------------------------------------
%  Flowcheck of the directed graph
%-------------------------------------------------------------------

%Goes through each combination of vertices. If they are not connected
%compute the amount of vertex distinct path from s to t and set
%MinConnectivity as the minimum of it
for i=1:AmountOfVertices
    for j=1:AmountOfVertices
        if AdjacencyMatrix(i,j)==0
            s=i;
            t=j;
            R=maxflow(G,s+AmountOfVertices,t);
            if R<MinConnectivity
                MinConnectivity=R;
                if MinConnectivity<3
                    %If the min connectivity is smaller than 3 we can
                    %stop
                    Result=false;
                    return
                end
            end
        end
    end
end

%otherwise the graph is 3-connected
Result=true;

end
