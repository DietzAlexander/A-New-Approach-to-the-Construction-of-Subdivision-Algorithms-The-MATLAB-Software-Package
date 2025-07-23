function [FaceMatrix] = computeFaceMatrix(AdjacencyMatrix,Inputcheck)

%COMPUTEFACEMATRIX Computes a face matrix out of a adjacency matrix belonging
%to a 3-connected planar graph by searching for induced cycles, which do not 
%separate the underlying graph. 
%
%Input: AdjacencyMatrix: n x n adjacency matrix
%       Inputcheck:      boolean. True, if an input check sould be done and
%                        false, if not
%Output: FaceMatrix: Matrix of faces. Each line represents a face with 
%indices of points in its entries. The size is amount of faces x maximal 
%face size 
%
%computeFaceMatrix(AdjacencyMatrix,Inputcheck) Computes a face matrix out of a
%adjacency matrix.


%Size of the adjacency matrix to check whether it is square
[rows,cols]=size(AdjacencyMatrix);
AmountOfVertices=rows;


%-------------------------------------------------------------------
%  Input Check
%-------------------------------------------------------------------

%Usually this function is part of other function in which the Inputcheck
%was already done. Therefore, it is optional.
if Inputcheck
    
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

    %Check for 3-vertex-connectivity
    if ~is3Connected(AdjacencyMatrix)
        error('Polygon is not 3-connected.');
    end

    %Check for planarity
    if ~PlanarityTest(AdjacencyMatrix)
        error('Polygon is not planar.');
    end


end




%-------------------------------------------------------------------
%  Computation of the faces
%-------------------------------------------------------------------

%Graph structure of the adjacency matrix
G=graph(AdjacencyMatrix);
GStart=G;

%Euler characteristic for planar graphs
AmountOfFaces=2-length(AdjacencyMatrix)+nnz(AdjacencyMatrix)/2;

%running variable for the while loop
FaceSize=3;

FaceMatrix=zeros(AmountOfFaces,height(G.Nodes));
EdgeBasisL=zeros(AmountOfFaces,height(G.Nodes));
EdgeBasisR=zeros(AmountOfFaces,height(G.Nodes));

EdgeList=G.Edges.EndNodes;
EdgeSeen=zeros(length(EdgeList),1);

BasisCounter=1;
while BasisCounter<=AmountOfFaces && FaceSize <= AmountOfVertices
    
    %The amount of faces found right now
    FoundFacesStartIteration=BasisCounter-1;
    
    %computes all cycles with a length of FaceSize
    [cycles,edgecycles]=allcycles(G,'MaxCycleLength',FaceSize,'MinCycleLength',FaceSize);

    %Goes trough all new cycles
    for i=1:height(edgecycles)
        currentEdgeCycle=cell2mat(edgecycles(i));
        currentCycle=cell2mat(cycles(i));
        

        %Checks if the cycle is an induced cycle
        InducedCheck=true;

        %Goes through all combinations of nodes in the cycle and check, if
        %there is a connection for 2 non neighbored nodes
        for k=1:length(currentCycle)
            for l=k+1:length(currentCycle)
                if k==1 && l==length(currentCycle)

                else
                    if l~=k+1
                        if AdjacencyMatrix(currentCycle(k),currentCycle(l))==1
                            InducedCheck=false;
                        end
                    end
                end
            end
        end

        %Graph without the cycle
        GTest=GStart;
        GTest=rmnode(GTest,currentCycle);

        %If the cycle is induced and the graph without the cycle is still
        %connected, it is a face and we add it to the basis list
        
        ConComp=conncomp(GTest);

        

        if InducedCheck && max(ConComp) == 1 
            FaceMatrix(BasisCounter,1:length(currentCycle))=currentCycle;
            EdgeListLocal=G.Edges.EndNodes(currentEdgeCycle,:);
            EdgeBasisL(BasisCounter,1:length(currentCycle))=min(EdgeListLocal,[],2);
            EdgeBasisR(BasisCounter,1:length(currentCycle))=max(EdgeListLocal,[],2);
            BasisCounter=BasisCounter+1;
        end

    end


    %In a 3-connected planar graphs an edge can be used in at most 2
    %faces. If an edge is used in 2 faces, it is removed from the graph to
    %speed the algorithm up (allcycles will run much faster)

    removeEdges=zeros(length(EdgeList),2);
    EdgeCounter=1;
    for i=FoundFacesStartIteration+1:BasisCounter-1
        EdgeListI=[EdgeBasisL(i,1:nnz(EdgeBasisL(i,:)));EdgeBasisR(i,1:nnz(EdgeBasisR(i,:)))]';


        [~,Position]=ismember(EdgeListI,EdgeList,'rows');
        NewRemovingEdges=EdgeListI(EdgeSeen(Position)==1,:);
        [AmountOfRemovingEdges,~]=size(NewRemovingEdges);
        removeEdges(EdgeCounter:EdgeCounter+AmountOfRemovingEdges-1,:)=NewRemovingEdges;
        EdgeSeen(Position)=1;

        EdgeCounter=EdgeCounter+AmountOfRemovingEdges;
    end
    if EdgeCounter>1
        G=rmedge(G,removeEdges(1:EdgeCounter-1,1),removeEdges(1:EdgeCounter-1,2));
    end
    FaceSize=FaceSize+1;

end

%-------------------------------------------------------------------
%  Tests if the results are valid
%-------------------------------------------------------------------

TestAdjacencyMatrix=FaceToAdjacencyMatrix(FaceMatrix);
if AdjacencyMatrix~=TestAdjacencyMatrix
    error('The computed face matrix differ from the original adjacency matrix. Was the graph not planar?')
end

%-------------------------------------------------------------------
%  Shrink the face matrix to a propper size
%-------------------------------------------------------------------

[row,~]=size(FaceMatrix);
MaxFace=0;
for i=1:row
    N=nnz(FaceMatrix(i,:));
    if N>MaxFace
        MaxFace=N;
    end
end

FaceMatrix=FaceMatrix(:,1:MaxFace);


end
