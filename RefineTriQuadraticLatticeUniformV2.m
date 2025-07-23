function [ControlpointsNew,LatticeNew,LatticeVolume,LatticeFaces,LatticeEdges,LatticeVertices,IsInnerVertex,ValenceInnerVertex] = RefineTriQuadraticLatticeUniformV2(Controlpoints,Lattice)
%REFINETRIQUADRATICLATTICEUNIFORMV2 refines a given lattice of control points 
% Controlpoints and volumes Lattice into a finer lattice using the 
% variant V2 of the tri quadratic subdivison algorithm.
%
%Input: Controlpoints in R^(n x 3)
%       Adjacency matrices of volumes Lattice in {0,1}^(n*n x AmountOfVolumes)
%       every column is an adjacency matrix with heigh and width of the
%       amount of control points. It can be restored with
%       reshape(P(:,i),[n,n])
%
%Output:    New controlpoints ControlpointsNew in the same format as the input
%           New volumes LatticeNew in the same format as the input
%            
%           Optional there is a separate output of the new volumes in the
%           categories PVolume, PFaces, PEdges, PVertices if needed (e.g.
%           for separate plotting). Also which old control point is an
%           inner vertex and its valence can be computed.
%


%-----------------------------------------------------------------
%       Input Check
%-----------------------------------------------------------------

%The amount of vertices and the dimension of the vertices
[AmountOfVertices,Dimension]=size(Controlpoints);


%The amount of volumes and the dimension of the adjacency matrices of the
%volumes
[AdjacencySize,AmountOfVolumes]=size(Lattice);

%The adjacency matrices needs to have a size in the amount of vertices
if sqrt(AdjacencySize)~=AmountOfVertices 
    error('Adjacency matrix has to have the size of vertices.')
end

%-----------------------------------------------------------------
%       Computation of general information for the algorithm
%-----------------------------------------------------------------

%Determines if a vetex is active on a volume. If the entry (i,j)==1 then
%the point i is active in the volume j
ActiveVertexInformation=zeros(AmountOfVertices,AmountOfVolumes);

%The amount of vertices per volume
VerticesPerVolume=zeros(AmountOfVolumes,1);

%The amount of edges per volume
EdgesPerVolume=zeros(AmountOfVolumes,1);

%The amount of faces per volume
FacesPerVolume=zeros(AmountOfVolumes,1);

%Goes through all volumes and compute the information
for i=1:AmountOfVolumes

    %The local adjacency matrix
    PLocal=reshape(Lattice(:,i),[AmountOfVertices,AmountOfVertices]);

    %If there is an entry >0 in a line of the adjacency matrix, the
    %correspondig point is active in the volume
    ActiveVertexInformation(:,i)=sum(PLocal,2)>0;

    %The sum of active points is the amount of points
    VerticesPerVolume(i)=sum(ActiveVertexInformation(:,i));

    % For each edge there are 2 entries in the adjacency matrix. So the sum
    % of all entries divides by 2 is the amount of edges
    EdgesPerVolume(i)=nnz(PLocal)/2;

    %Use the euler characteristic to get the amount of faces per volume
    FacesPerVolume(i)=2+EdgesPerVolume(i)-VerticesPerVolume(i);
end

%Summing up all vertices in all volumes is the amount of new points (as
%every volume is mapped to a smaller volume
AmountOfNetPoints=sum(sum(ActiveVertexInformation));
ControlpointsNew=zeros(AmountOfNetPoints,Dimension);

%For each volume there is a new volume created. This variable contains the
%information, which vertex is the first vertex in the list of all vertices
%in the corresponing volume (to get an easier identification)
FirstVertexOfVolumeInVNew=zeros(AmountOfVolumes,1);

%The point from which the new point was refined
ParentPoint=zeros(AmountOfNetPoints,1);



%-----------------------------------------------------------------
%       Points & Volumes out of Volumes
%-----------------------------------------------------------------

%The indices and values of the new adjacency matrix of the new volumes out
%of volumes. As every volume is refinded to a smaller one of same size, the
%adjacency matrices are identical up to the numbers ob rows and columns.
%As for sparse matrices the indexing is much more faster, we do it in this
%way
Prelocation=nnz(Lattice);

IndicesRow=zeros(Prelocation,1);
IndicesCol=zeros(Prelocation,1);


%Identification value for point numbers
PointList=1:AmountOfNetPoints;

%Counter for the indices
IndexCounter=1;

%Counter for the new points
PointCounter=1;

%Goes through all volumes. Creates the new points and the new adjacency
%matrices for the volumes out of volumes
for i=1:AmountOfVolumes

    %Produces the  small or condensated adjacency matrix for a volume

    %First the current matrix is put in the right shape
    PLocal=reshape(Lattice(:,i),[AmountOfVertices,AmountOfVertices]);

    %Second get just the active lines
    PSmall=PLocal(logical(ActiveVertexInformation(:,i)),logical(ActiveVertexInformation(:,i)));

    %Computes the subdivision matrix to refine vertices
    SubdivisionMatrixSmall=computeTriQuadraticSubdivisionMatrixV2(PSmall);

    %Refine the vertices
    VNewLocal=SubdivisionMatrixSmall*Controlpoints(logical(ActiveVertexInformation(:,i)),:);

    %Set the new vertices to the vector of all vertices
    ControlpointsNew(PointCounter:PointCounter+length(VNewLocal)-1,:)=VNewLocal;

    %Set the value of the first vertex of the volume
    FirstVertexOfVolumeInVNew(i)=PointCounter;

    %Set the parent information for each new point
    ParentPoint(PointCounter:PointCounter+length(VNewLocal)-1)=PointList(logical(ActiveVertexInformation(:,i)));

    %Compute the new adjacency matrix for each volume out of volume

    %Get the indices of the non zero entries of the old one
    [rowLocal,colLocal]=find(PSmall);

    %Get the amount of entries for each volume
    AmountOfNewEntries=length(rowLocal);

    %Set the row index by a computation out of the col and the row vector
    %(the final matrix has just one column for the whole matrix so the index
    %computation is a little more difficult)
    IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=(PointCounter+colLocal-1-1)*AmountOfNetPoints+PointCounter+rowLocal-1;

    %The column index is just the i-th column)
    IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=i;

    %raise the two new counter
    PointCounter=PointCounter+length(VNewLocal);
    IndexCounter=IndexCounter+AmountOfNewEntries;
end

%The variable for the new volume information (adjacency matrices of each volume)
LatticeVolume=sparse(IndicesRow,IndicesCol,ones(length(IndicesRow),1),AmountOfNetPoints*AmountOfNetPoints,AmountOfVolumes);

%-----------------------------------------------------------------
%       Preparation for Face Identification
%-----------------------------------------------------------------

%-----------------------------------------------------------------
%In a first step we compute all adjacency matrices of all faces in all old
%volumes

%The maximum amount of entries in the face adjacency matrix. If every face
%consists out of all vertices
Prelocation=2*max(VerticesPerVolume)*sum(FacesPerVolume);

%Index for the following sparse matrix
IndicesRow=zeros(Prelocation,1);
IndicesCol=zeros(Prelocation,1);

%Index for the following sparse matrix
IndicesRowRefined=zeros(Prelocation,1);
IndicesColRefined=zeros(Prelocation,1);


%Counter for the indices
IndexCounter=1;
FaceCounter=1;

%Goes through all volumes and compute the adjacency matrix of the
%corresponding faces
for i=1:AmountOfVolumes

    %Compute the local adjacency matrix
    PLocal=reshape(Lattice(:,i),[AmountOfVertices,AmountOfVertices]);

    %Delete all empty lines 
    PSmall=PLocal(logical(ActiveVertexInformation(:,i)),logical(ActiveVertexInformation(:,i)));

    %The indices of the current active points.
    ActivePointsLocal=find(ActiveVertexInformation(:,i));

    %Get the current active points of this volume
    AmountOfActivePointsLocal=sum(ActiveVertexInformation(:,i));

    %Compute the faces of the local volume
    FacesLocal=computeFaceMatrix(PSmall,false);

    %Get the size of the largest face of the volume 
    [~,maxCol]=size(FacesLocal);

    %Goes through all faces of the volume and compute the corresponding
    %face matrix
    for j=1:FacesPerVolume(i)

        %compute the local version of th adjacency matrix
        FaceAdjacencyMatrixSmall=zeros(AmountOfActivePointsLocal,AmountOfActivePointsLocal);
        
        %Write the entries of the local adjacency matrix of the face
        for k=1:maxCol
            if FacesLocal(j,k)>0
                if k==maxCol || FacesLocal(j,k+1)==0
                    FaceAdjacencyMatrixSmall(FacesLocal(j,k),FacesLocal(j,1))=1;
                    FaceAdjacencyMatrixSmall(FacesLocal(j,1),FacesLocal(j,k))=1;
                else
                    FaceAdjacencyMatrixSmall(FacesLocal(j,k),FacesLocal(j,k+1))=1;
                    FaceAdjacencyMatrixSmall(FacesLocal(j,k+1),FacesLocal(j,k))=1;
                end
            end
        end


        %Get the indices of the non zero entries of the old one
        [rowLocal,colLocal]=find(FaceAdjacencyMatrixSmall);

        %Get the amount of entries for each volume
        AmountOfNewEntries=length(rowLocal);

        %Set the row index by a computation out of the col and the row vector
        %(the final matrix has just one column for the whole matrix so the index
        %computation is a little more difficult)
        IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=(ActivePointsLocal(colLocal)-1)*AmountOfVertices+ActivePointsLocal(rowLocal);

        %The column index is just the face counter
        IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=FaceCounter;

        %Set the row index by a computation out of the col and the row vector
        %(the final matrix has just one column for the whole matrix so the index
        %computation is a little more difficult)
        IndicesRowRefined(IndexCounter:IndexCounter+AmountOfNewEntries-1)=(FirstVertexOfVolumeInVNew(i)+colLocal-1-1)*AmountOfNetPoints+FirstVertexOfVolumeInVNew(i)+rowLocal-1;
    
        %The column index is just the face counter
        IndicesColRefined(IndexCounter:IndexCounter+AmountOfNewEntries-1)=FaceCounter;

        %Raising the counter
        FaceCounter=FaceCounter+1;
        IndexCounter=IndexCounter+AmountOfNewEntries;

    end
end

%Computes the matrix of all adjacency matrices of all faces. Here it is
%transposed, so every line is an adjacency matrix (preperation for the
%unique command)
FaceAdjacencyMatrixT=sparse(IndicesCol(1:IndexCounter-1),IndicesRow(1:IndexCounter-1),ones(IndexCounter-1,1),FaceCounter-1,AmountOfVertices*AmountOfVertices);

%Computes the matrix of all adjacency matrices of all refined faces. Here it is
%transposed, so every line is an adjacency matrix (preperation for the
%unique command)
FaceAdjacencyMatrixRefinedT=sparse(IndicesColRefined(1:IndexCounter-1),IndicesRowRefined(1:IndexCounter-1),ones(IndexCounter-1,1),FaceCounter-1,AmountOfNetPoints*AmountOfNetPoints);

%Get a list of all indizes of the non zero entries of the refinded faces. 
%They are needed for a fast setting if indizes for the following sparse
%matrices
[RowFaceAdjacencyMatrixUniqueRefindedT,ColFaceAdjacencyMatrixUniqueRefinedT]=find(FaceAdjacencyMatrixRefinedT);


%The amount of comuted faces
[AmountOfComputedFaces,~]=size(FaceAdjacencyMatrixT);

%-----------------------------------------------------------------
%We compute now an addition to the face matrix, that every line has the
%same amount of ones

%The amount of non zero entries of every face
NNZPerFace=sum(FaceAdjacencyMatrixT,2);


%The maximal amount of non zero entries of every face
NMax=max(NNZPerFace);

%The maximal amount of non zero entries of every face
Nmin=min(NNZPerFace);


%The maximum amount of entries in the face adjacency matrix. If every face
%consists out of all vertices
Prelocation=AmountOfComputedFaces*NMax-sum(NNZPerFace);

%Index for the following sparse matrix
IndicesRow=zeros(Prelocation,1);
IndicesCol=zeros(Prelocation,1);


IndexCounter=1;

%Goes through all computed faces and adds ones to its adjacency matrix. The
%reason for that is, that we can easyly get and reduce the nnz information
%of each line, if each line has the same amount of nnz.
for i=1:AmountOfComputedFaces

    %The amount of new entries (maximal face nnz - current face nnz)
    AmountOfNewEntries=NMax-NNZPerFace(i);

    %Adds the indices of entries that become 1.
    IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=i*ones(NMax-NNZPerFace(i),1);
    IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=1:(NMax-NNZPerFace(i));

    %raise the counter
    IndexCounter=IndexCounter+AmountOfNewEntries;
    
end

%Creates the matrix of additional ones
AdditionMatrix=sparse(IndicesRow(1:IndexCounter-1),IndicesCol(1:IndexCounter-1),ones(IndexCounter-1,1),AmountOfComputedFaces,NMax-Nmin);

%Add the ones to the face adjacency matrix
FaceAdjacencyMatrixAddT=sparse([FaceAdjacencyMatrixT,AdditionMatrix]);

%-----------------------------------------------------------------
%We now unique the face adjacency matrix to get later the information of
%the inner faces

%Reduces the adjacency matrix to just information about the ones. That is a
%good coice because the unique command do not have to deal with the
%adjacency matrix but just with this reduced version and is therefore much
%faster
[row,~]=find(FaceAdjacencyMatrixAddT');
ReducteInformation=reshape(row,[NMax,AmountOfComputedFaces])';

%Uniques the faces
[~,FaceList,FaceIdentifyer]=unique(ReducteInformation,'rows','stable');

%Get the uniqued face adjacency matrices
FaceAdjacencyMatrixUniqueT=FaceAdjacencyMatrixT(FaceList,:);

%-----------------------------------------------------------------
%We now compute the active points on each refined faces to compute later
%the adjacency matrices for volumes out of faces, edges and vertices

%The amount of refinded faces
[AmountOfRefinedFaces,~]=size(FaceAdjacencyMatrixRefinedT);

%The variable containing the active points per refined face
ActivePointsInRefinedFace=zeros(AmountOfRefinedFaces,NMax/2);

%Goes through  all refined faces
for i=1:AmountOfRefinedFaces

    %Get the column values for every non zero entry of the adjacency
    %matrix. As the adjacency matrix is coded in lines the values are
    %between 1 and AmountOfNetPoints*AmountOfNetPoints
    Line=ColFaceAdjacencyMatrixUniqueRefinedT(RowFaceAdjacencyMatrixUniqueRefindedT==i);

    %Reduce the information to the columns of the square matrix version
    Line=mod(Line,AmountOfNetPoints);

    %Uniques the columns
    ActivePoints=unique(Line);

    %Set the value for the last column (as modulo sets it to 0)
    ActivePoints(ActivePoints==0)=AmountOfNetPoints;

    %Put the value into the matrix
    ActivePointsInRefinedFace(i,1:length(ActivePoints))=ActivePoints;
end

%-----------------------------------------------------------------
%       Volumes out of faces
%-----------------------------------------------------------------

%The amount of inner faces
AmountOfInnerFaces=0;

%Determines if a face is an inner face or not
IsInnerFace=zeros(length(FaceList),1);

%The amount of reduced (uniqued) faces
AmountOfReducedFaces=length(FaceList);

%Prelocation for the face matrix If a face has NMax non zero values, there 
% are 3 times nnz values in the corresponding volume (top, bottom, lines between) 
Prelocation=AmountOfReducedFaces*NMax*3;

%Index for the following sparse matrix. 
IndicesRow=zeros(Prelocation,1);
IndicesCol=zeros(Prelocation,1);

%Index counter
IndexCounter=1;

%Goes through all reduced (uniqued) faces
for i=1:AmountOfReducedFaces

    %Determines which one of the refinded faces correspond with this face
    RefinedFaceIdentifyer=find(FaceIdentifyer==i);

    %If there is just one refined face, it is an outer face
    if isscalar(RefinedFaceIdentifyer)
        
    %If there are more than 2 refined faces, it is not well defined
    elseif length(RefinedFaceIdentifyer)>2
        error('At most 2 volumes has to have a common face.')
    
    %Else it is an inner face
    else

        %Set the inner face value to 1
        IsInnerFace(i)=1;

        %Raise the amount of inner faces
        AmountOfInnerFaces=AmountOfInnerFaces+1;

        %-----------------------------------------------------------------
        %Top and bottom face 

        %Get the indices of the top and bottom refined face (which are both
        %faces of the new volume) and adds it to the adjacency matrix of
        %the new volume
        rowLocal=[ColFaceAdjacencyMatrixUniqueRefinedT(RowFaceAdjacencyMatrixUniqueRefindedT==RefinedFaceIdentifyer(1));...
                    ColFaceAdjacencyMatrixUniqueRefinedT(RowFaceAdjacencyMatrixUniqueRefindedT==RefinedFaceIdentifyer(2))];

        %Get the amount of new entries
        AmountOfNewEntries=length(rowLocal);

        %Set the row index 
        IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=rowLocal;

        %The column index is the new volume (so the actual counter of inner
        %faces
        IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=AmountOfInnerFaces;

        %raise the index counter
        IndexCounter=IndexCounter+AmountOfNewEntries;

        %-----------------------------------------------------------------
        %Side faces

        %Get the active points of the top face
        AmountOfActivePoints1=nnz(ActivePointsInRefinedFace(RefinedFaceIdentifyer(1),:));
        ActivePoints1=ActivePointsInRefinedFace(RefinedFaceIdentifyer(1),1:AmountOfActivePoints1)';

        %Get the active points of the bottom face
        AmountOfActivePoints2=nnz(ActivePointsInRefinedFace(RefinedFaceIdentifyer(2),:));
        ActivePoints2=ActivePointsInRefinedFace(RefinedFaceIdentifyer(2),1:AmountOfActivePoints2)';

        %Put both together
        ActiveLocalPoints=[ActivePoints1;ActivePoints2];

        %and get the parent information
        LocalParents=ParentPoint(ActiveLocalPoints);

        %Goes through all vertices of one face and find the corresponding
        %other point
        for j=1:length(LocalParents)/2

            %Connect both with a line
            Connection=find(LocalParents==LocalParents(j));
            rowLocal=[ActiveLocalPoints(Connection(1));ActiveLocalPoints(Connection(2))];
            colLocal=[ActiveLocalPoints(Connection(2));ActiveLocalPoints(Connection(1))];

            %and set the adjacency values

            %Get the amount of entries for each volume
            AmountOfNewEntries=length(rowLocal);
    
            %Set the row index by a computation out of the col and the row vector
            %(the final matrix has just one column for the whole matrix so the index
            %computation is a little more difficult)
            IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=(colLocal-1)*AmountOfNetPoints+rowLocal;
    
            %The column index is the new volume (so the actual counter of inner
            %faces
            IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=AmountOfInnerFaces;
            
            %raise the index counter
            IndexCounter=IndexCounter+AmountOfNewEntries;
        end
   
    end

end

%The variable for the new volume information (adjacency matrices of each volume out of a face)
LatticeFaces=sparse(IndicesRow(1:IndexCounter-1),IndicesCol(1:IndexCounter-1),ones(IndexCounter-1,1),AmountOfNetPoints*AmountOfNetPoints,AmountOfInnerFaces);



%-----------------------------------------------------------------
%       Volume out of an Edge
%-----------------------------------------------------------------

%-----------------------------------------------------------------
%Starts by getting all edges of the (not refined) lattice)

%Creates an adjacency matrix for the whole graph
WholeGraphAdjacency=sum(Lattice,2);
WholeGraphAdjacency(WholeGraphAdjacency>1)=1;

%Creates the whole graph
WholeGraph=graph(reshape(WholeGraphAdjacency,[AmountOfVertices,AmountOfVertices]));

%Get all edges of the whole graph
WholeEdges=WholeGraph.Edges.EndNodes;

%Information if an edge is an inner edge
IsInnerEdge=zeros(length(WholeEdges),1);

%Valence of the corresponding edge
ValenceEdge=zeros(length(WholeEdges),1);

%Prelocation for the edge matrix. For every edge on every inner face (there
%are doubles of course) there are 8 entries created (2*4 as every edge is
%double represented)
Prelocation=nnz(FaceAdjacencyMatrixUniqueT)/2*8;

%Index for the following sparse matrix
IndicesRow=zeros(Prelocation,1);
IndicesCol=zeros(Prelocation,1);
Entries=zeros(Prelocation,1);

%Index counter
IndexCounter=1;

%The amount of inner edges
AmountOfInnerEdges=0;

%Goes through all edges
for i=1:length(WholeEdges)

    %Set the current edge
    CurrentEdge=WholeEdges(i,:);

    %Get the entry if the edge (in adjacency matrix). There are 2 but one of
    %them is enough
    IdentifyerInAdjacencyMatrix=(CurrentEdge(1)-1)*AmountOfVertices+CurrentEdge(2);

    %Get the amount of faces on this edge
    FaceOnCurrentEdge=find(FaceAdjacencyMatrixUniqueT(:,IdentifyerInAdjacencyMatrix));

    %Amount if faces on the current edge
    ValenceEdge(i)=length(FaceOnCurrentEdge);

    %Checks if the edge is an inner edge (which is if every face having
    %this edge is an inner face
    if sum(IsInnerFace(FaceOnCurrentEdge))==length(FaceOnCurrentEdge)

        %Set the inner edge information
        IsInnerEdge(i)=1;

        %Raise the amount of inner edge
        AmountOfInnerEdges=AmountOfInnerEdges+1;

        %Goes through all inner faces and produce the side face of the
        %prism
        for j=1:ValenceEdge(i)

            %Get the current face
            RefinedFaceIdentifyer=find(FaceIdentifyer==FaceOnCurrentEdge(j));
            
            %Get the current edge of the 2 refined volumes with the inner
            %face
            AmountOfActivePoints1=nnz(ActivePointsInRefinedFace(RefinedFaceIdentifyer(1),:));
            ActivePoints1=ActivePointsInRefinedFace(RefinedFaceIdentifyer(1),1:AmountOfActivePoints1)';
    
            AmountOfActivePoints2=nnz(ActivePointsInRefinedFace(RefinedFaceIdentifyer(2),:));
            ActivePoints2=ActivePointsInRefinedFace(RefinedFaceIdentifyer(2),1:AmountOfActivePoints2)';

            %Get the indices of the two points of the first volume
            ActivePoint1A=ActivePoints1(ParentPoint(ActivePoints1)==CurrentEdge(1));
            ActivePoint1B=ActivePoints1(ParentPoint(ActivePoints1)==CurrentEdge(2));

            %Get the indices of the two points of the second volume
            ActivePoint2A=ActivePoints2(ParentPoint(ActivePoints2)==CurrentEdge(1));
            ActivePoint2B=ActivePoints2(ParentPoint(ActivePoints2)==CurrentEdge(2));

            %Creates the 8 entries:
            %Edge in first volume both directions
            %Edge in second volume both directions
            %Connection of first point both directions
            %Connection of second point both directions

            LocalRow=[ActivePoint1A;...
                    ActivePoint1B;...
                    ActivePoint2A;...
                    ActivePoint2B;...
                    ActivePoint1A;...
                    ActivePoint2A;...
                    ActivePoint1B;...
                    ActivePoint2B];

            LocalCol=[ActivePoint1B;...
                    ActivePoint1A;...
                    ActivePoint2B;...
                    ActivePoint2A;...
                    ActivePoint2A;...
                    ActivePoint1A;...
                    ActivePoint2B;...
                    ActivePoint1B];
            
            %The first 4 enrties are just 1/2 as the edge in each volume is
            %visited twice
            LocalEntries=[1/2;...
                        1/2;...
                        1/2;...
                        1/2;...
                        1;...
                        1;...
                        1;...
                        1];

            %Get the amount of entries for each new volume
            AmountOfNewEntries=length(LocalRow);
    
            %Set the row index by a computation out of the col and the row vector
            %(the end matrix has just one column for the whole matrix so the index
            %computation is a little more difficult)
            IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=(LocalCol-1)*AmountOfNetPoints+LocalRow;
    
            %The column index is just the edge counter
            IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=AmountOfInnerEdges;

            %Set the entries
            Entries(IndexCounter:IndexCounter+AmountOfNewEntries-1)=LocalEntries;

            %Raise the index counter
            IndexCounter=IndexCounter+AmountOfNewEntries;
        end

        

    end
end

%The variable for the new volume information (adjacency matrices of each volume out of an edge)
LatticeEdges=sparse(IndicesRow(1:IndexCounter-1),IndicesCol(1:IndexCounter-1),Entries(1:IndexCounter-1),AmountOfNetPoints*AmountOfNetPoints,AmountOfInnerEdges);



%-----------------------------------------------------------------
%       Volume out of a point
%-----------------------------------------------------------------


%Information if the vertex is an inner vertex
IsInnerVertex=zeros(AmountOfVertices,1);

%Valence of the inner vertex (Amount of faces having this vertex in common)
ValenceInnerVertex=zeros(AmountOfVertices,1);

%The amount of edges in every face is also the amount of vertices in every
%face. For each face a vertex is in there are set one edge, so 2
%information in the adjacency matrix
Prelocation=nnz(FaceAdjacencyMatrixUniqueT)/2*2;

%Index for the following sparse matrix
IndicesRow=zeros(Prelocation,1);
IndicesCol=zeros(Prelocation,1);

%The index counter
IndexCounter=1;

%The amount of inner vertices
AmountOfInnerVertices=0;

%Goes through all vertices 
for i=1:AmountOfVertices

    %Find all faces having this vertex on common
    [FaceOnCurrentVertex,~]=find(FaceAdjacencyMatrixUniqueT(:,AmountOfVertices*(i-1)+1:AmountOfVertices*(i)));

    %reduces the information to the face numbers
    FaceOnCurrentVertex=unique(FaceOnCurrentVertex);

    %Set the valence of the vertex as amount of faces who have this vertex
    %in common
    ValenceInnerVertex(i)=length(FaceOnCurrentVertex);

    %checks if the vertex is an inner vertex (that is if every face having
    %this vertex is an inner face
    if sum(IsInnerFace(FaceOnCurrentVertex))==length(FaceOnCurrentVertex)

        %Set the information, that the vertex is an inner vertex
        IsInnerVertex(i)=1;
            
        %Raise the amount of inner vertices
        AmountOfInnerVertices=AmountOfInnerVertices+1;

        %Goes through all faces having this vertex in common
        for j=1:ValenceInnerVertex(i)

            %Get the refined volumes having this face in common
            RefinedFaceIdentifyer=find(FaceIdentifyer==FaceOnCurrentVertex(j));
            
            %Get the active points of each refinded face
            AmountOfActivePoints1=nnz(ActivePointsInRefinedFace(RefinedFaceIdentifyer(1),:));
            ActivePoints1=ActivePointsInRefinedFace(RefinedFaceIdentifyer(1),1:AmountOfActivePoints1)';

            AmountOfActivePoints2=nnz(ActivePointsInRefinedFace(RefinedFaceIdentifyer(2),:));
            ActivePoints2=ActivePointsInRefinedFace(RefinedFaceIdentifyer(2),1:AmountOfActivePoints2)';

            %Get the point in the refined face, that belongs to the current
            %vertex
            ActivePoint1=ActivePoints1(ParentPoint(ActivePoints1)==i);


            ActivePoint2=ActivePoints2(ParentPoint(ActivePoints2)==i);

            %Set the indices for the adjacency matrix
            LocalRow=[ActivePoint1;...
                    ActivePoint2];

            LocalCol=[ActivePoint2;...
                    ActivePoint1];


            %Get the amount of entries for each volume
            AmountOfNewEntries=length(LocalRow);
    
            %Set the row index by a computation out of the col and the row vector
            %(the final matrix has just one column for the whole matrix so the index
            %computation is a little more difficult)
            IndicesRow(IndexCounter:IndexCounter+AmountOfNewEntries-1)=(LocalCol-1)*AmountOfNetPoints+LocalRow;
    
            %The column index is just the current volume belonging to the
            %vertex
            IndicesCol(IndexCounter:IndexCounter+AmountOfNewEntries-1)=AmountOfInnerVertices;

            %Raise the index counter
            IndexCounter=IndexCounter+AmountOfNewEntries;
        end
    end
end

%The variable for the new volume information (adjacency matrices of each volume out of a vertex)
LatticeVertices=sparse(IndicesRow(1:IndexCounter-1),IndicesCol(1:IndexCounter-1),ones(IndexCounter-1,1),AmountOfNetPoints*AmountOfNetPoints,AmountOfInnerVertices);


%Put all four types of new volumes together
LatticeNew=sparse([LatticeVolume,LatticeFaces,LatticeEdges,LatticeVertices]);

%If there occur double counts at the indices, set the values to 1 to get
%adjacency matrices
LatticeNew(LatticeNew>1)=1;


end
