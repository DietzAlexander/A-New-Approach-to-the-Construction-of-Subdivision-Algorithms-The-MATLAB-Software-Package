function [Result] = PlanarityTest(AdjacencyMatrix)
%PLANARITYTEST Checks, if the graph related to the AdjacencyMatrix is 
%planar or not. Returns true if planar and false, if not.
%
%Input: n x n adjacency matrix
%Output: Boolean true or false
%
%PlanarityTest(AdjacencyMatrix) checks if the graph to AdjacencyMatrix is
%planar.


%During the algorithm this value gets false, if the graph is non planar and
%the algorithm stops. If it stays true and there is no return up to the end
%of the algorithm, the Result stays true and is given back.
Result=true;

%-------------------------------------------------------------------
%  Input Check
%-------------------------------------------------------------------

%Size of the adjacency matrix to check whether it is square
[rows,cols]=size(AdjacencyMatrix);
AmountOfVertices=rows;

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
%  Necesarry condition for planarity
%-------------------------------------------------------------------

%Check necesarry condition for planarity
if nnz(AdjacencyMatrix)/2>3*length(AdjacencyMatrix)-6
    Result=false;
    return
end


%-------------------------------------------------------------------
%  Creates variables for the DFS tree
%-------------------------------------------------------------------

%Graph representation of the adjacency matrix
G=graph(AdjacencyMatrix);

%List of the edges (during the algorithm eges have a direction, which is
%coded in this list (every edge goes from EdgeList(1) to EdgeList(2).
EdgeList=G.Edges.EndNodes;

%The parent edge of a vertex in the DFS tree
ParentEdge=zeros(G.numnodes,1);

%The height of a vertex in the DFS tree (0 is the root)
Height=ones(G.numnodes,1)*G.numnodes+1;

%A value of an edge. Describes the lowest vertex, which can be reached with
%this edge or any path starting with this edge in the directed graph of the
%DFS tree and its back edges. The value is Edge(2) for each back edge.
LowPoint=zeros(G.numedges,1);

%Value for the first DFS algorithm which determines, if an edge was seen or
%used resp.
EdgeUsed=zeros(G.numedges,1);

%Describes if an edge is a tree edge or a back edge
IsBackEdge=zeros(G.numedges,1);

Height(1)=0;

%-------------------------------------------------------------------
%  Creates a DFS tree and the LowPoint value of each edge.
%-------------------------------------------------------------------

%Creates a DFS tree and the Lowpoint value of each edge.
[EdgeList,ParentEdge,~,LowPoint,~,IsBackEdge]=DFS1(1,EdgeList,ParentEdge,Height,LowPoint,EdgeUsed,IsBackEdge,G);

%-------------------------------------------------------------------
%  Preprocessing for the creation of the constraint matrix
%-------------------------------------------------------------------

%The constraint matrix
ConstraintMatrix=zeros(G.numedges);

%Creates the fundamental cycle for each back edge (A cycle which consists
%of the back edge and just tree edges)
CycleList=zeros(G.numedges);
CycleListLength=zeros(G.numedges,1);
for i=1:G.numedges
    %Just look at the edge, if it is a back edge 
    if IsBackEdge(i)==1
        %Compute the fundamental cycle of the edge i
        CurrentEdge=i;
        CycleList(i,1)=i;
        counter=2;
        Target=EdgeList(CurrentEdge,1);
        while Target ~= EdgeList(i,2)
            CurrentEdge=ParentEdge(Target);
            CycleList(i,counter)=CurrentEdge;
            Target=EdgeList(CurrentEdge,1);
            counter=counter+1;
        end
        CycleListLength(i)=nnz(CycleList(i,:));
    end
end

%Creates all edges, that are connected to a node
NodeEdgeList=zeros(G.numnodes,G.numedges);
NodeEdgeListLength=zeros(G.numnodes,1);

for i=1:G.numnodes
    [row,~]=find(EdgeList==i);
    NodeEdgeList(i,1:length(row))=row;
    NodeEdgeListLength(i)=length(row);
end


%-------------------------------------------------------------------
%  Creation of the constraint matrix
%-------------------------------------------------------------------

%Compares each pair of back edges if they have same and opposite
%constraints
for i=1:G.numedges
    %Just look at the edge, if it is a back edge 
    if IsBackEdge(i)==1

        %Compute the fundamental cycle of the edge i
        CycleListI=CycleList(i,1:CycleListLength(i));

        %Go through all edges greater than the current (avoid double check)
        for j=i+1:G.numedges
            %Just look at the edge, if it is a back edge 
            if IsBackEdge(j)==1
                CycleListJ=CycleList(j,1:CycleListLength(j));

                %Look at the constrains if both edges has an overlapping
                %fundamental cycle
                if sum(ismember(CycleListI,CycleListJ))>0
                    
                    E1=CycleListI(find(ismember(CycleListI,CycleListJ),1)-1);
                    E2=CycleListJ(find(ismember(CycleListJ,CycleListI),1)-1);

                    %Check for a different constraint
                    if LowPoint(E2)<LowPoint(i) && LowPoint(E1)<LowPoint(j)
                        ConstraintMatrix(i,j)=-1;
                        ConstraintMatrix(j,i)=-1;
                    end

                    %Check for a self constraint
                    SameEdges=CycleListI(ismember(CycleListI,CycleListJ));
                    CompareValue=min(LowPoint(i),LowPoint(j));
                    for k=1:length(SameEdges)-1
                        vk=EdgeList(SameEdges(k),1);
                        row=NodeEdgeList(vk,1:NodeEdgeListLength(vk));
                        for l=1:length(row)
                            if row(l)~=SameEdges(k) && row(l)~=SameEdges(k+1) 
                                if LowPoint(row(l))<CompareValue
                                    if ConstraintMatrix(i,j)==-1
                                        Result=false;
                                        return
                                    else
                                        ConstraintMatrix(i,j)=1;
                                        ConstraintMatrix(j,i)=1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


%Deletes every edge, that has no constraint in the constraint matrix
i=1;
while i<=length(ConstraintMatrix)
    if nnz(ConstraintMatrix(i,:))==0 && nnz(ConstraintMatrix(:,i))==0
        ConstraintMatrix(i,:)=[];
        ConstraintMatrix(:,i)=[];
    else
        i=i+1;
    end
end

%-------------------------------------------------------------------
%  Check if the constraint graph is balanced
%-------------------------------------------------------------------

SignGraph=graph(ConstraintMatrix);
ConnectedComponent=conncomp(SignGraph);

%Check for every connected component
for i=1:max(ConnectedComponent)
    %Compute the local constraint matrix
    LocalSignGraph=graph(ConstraintMatrix(ConnectedComponent==i,ConnectedComponent==i));

    %Initialte values for the DFS
    NodeSeen=zeros(LocalSignGraph.numnodes,1);
    NodeSign=zeros(LocalSignGraph.numnodes,1);
    IsTreeEdge=zeros(LocalSignGraph.numedges,1);
    v=1;
    NodeSign(v)=1;
    NodeSeen(1)=1;
    EdgeList=LocalSignGraph.Edges.EndNodes;
    EdgeSign=LocalSignGraph.Edges.Weight;
    
    %Compute the DFS and give every node a + or - sign
    [~,NodeSign,IsTreeEdge,EdgeList,EdgeSign] = DFSSign(NodeSeen,NodeSign,IsTreeEdge,EdgeList,EdgeSign,v);

    %Goes trought all non tree edges
    for j=1:LocalSignGraph.numedges
        if IsTreeEdge(j)==0
            %Check, if the sign of the back edge is valid or not
            if NodeSign(EdgeList(j,1))*NodeSign(EdgeList(j,2))~=EdgeSign(j)
                Result=false;
                return
            end
        end
    end


end


end







%Does a DFS and gives every node a + or - sign
function [NodeSeen,NodeSign,IsTreeEdge,EdgeList,EdgeSign] = DFSSign(NodeSeen,NodeSign,IsTreeEdge,EdgeList,EdgeSign,v)
    [row,col]=find(EdgeList==v);
    for i=1:length(row)
        w=EdgeList(row(i),mod(col(i),2)+1);
        if NodeSeen(w)==0
            IsTreeEdge(row(i))=1;
            NodeSeen(w)=1;
            NodeSign(w)=NodeSign(v)*EdgeSign(row(i));
            [NodeSeen,NodeSign,IsTreeEdge,EdgeList,EdgeSign] = DFSSign(NodeSeen,NodeSign,IsTreeEdge,EdgeList,EdgeSign,w);
        end
    end

end


%Computes recursively the DFS tree of the graph and the lowpoint value of
%each edge.
function [EdgeList,ParentEdge,Height,LowPoint,EdgeUsed,IsBackEdge] = DFS1(v,EdgeList,ParentEdge,Height,LowPoint,EdgeUsed,IsBackEdge,G)
    %Get the parent edge, if possible
    if Height(v)~=0
        e=ParentEdge(v);
    else
        e=0;
    end
    
    %Search for edges not visited right now
    [row,~]=find(EdgeList==v);

    %Go throught all non visited edges
    while sum(EdgeUsed(row))<length(row)

        %Select one unused edge
        UnusedEdges=row((EdgeUsed(row)==0));
        CurrentEdge=UnusedEdges(1);
        EdgeUsed(CurrentEdge)=1;

        %Code the orientation of the edge and swap it if nececarry
        if EdgeList(CurrentEdge,1)~=v
            z=EdgeList(CurrentEdge,1);
            EdgeList(CurrentEdge,1)=EdgeList(CurrentEdge,2);
            EdgeList(CurrentEdge,2)=z;
        end
        
        LowPoint(CurrentEdge)=Height(v);
        
        w=EdgeList(CurrentEdge,2);

        %Edge is a tree edge
        if Height(w)==G.numnodes+1
            ParentEdge(w)=CurrentEdge;
            Height(w)=Height(v)+1;
            [EdgeList,ParentEdge,Height,LowPoint,EdgeUsed,IsBackEdge] = DFS1(w,EdgeList,ParentEdge,Height,LowPoint,EdgeUsed,IsBackEdge,G);    
            
        %Edge is a back edge
        else
            IsBackEdge(CurrentEdge)=1;
            LowPoint(CurrentEdge)=Height(w);
        end

        if e~=0
            if LowPoint(CurrentEdge)<LowPoint(e)
                LowPoint(e)=LowPoint(CurrentEdge);
            end
        end
    end
end

