function [ControlpointsNew,LatticeNew] = RefineTriCubicLatticeUniform(Controlpoints,Lattice)
% REFINETRUCUBICLATTICEUniform refines the given cubic lattice globally.
% Each initial element consisting out of cubes sharing an inner vertex is
% refined to a smaller one by using the cubic subdivision algorithm.
%
% Input:    Lattice:            The lattice. Each column stands for one cube.
% Output:   LatticeNew:         The new lattice
%           ControlpointsNew:   The new control points
%
% RefineTriCubicLatticeStructure(Lattice) refines the given cubic lattice globally. 

%The amount of vertices and the amount of cubes
[AmountOfVertices,AmountOfVolumes]=size(Lattice);
AmountOfVertices=sqrt(AmountOfVertices);



%The amount of faces
AmountOfFaces=sum((sum(Lattice)/2)/2+2-2);

%As all volumes are cubes, all faces are quads
MaxFaceSize=4;



%The maximum amount of faces per volume (complicated formula for 6)
MaxFacesPerVolume=max((sum(Lattice)/2)/2+2-2);

%The maximum amount of faces per volume (complicated formula for 8)
MaxVerticesPerVolume= max(sum(Lattice)/2/2+2);

%-------------------------------------------------------------------------
%       Compute the information of the faces of the cubes
%-------------------------------------------------------------------------

%Two variables to avaiod, that computeFaceMatrix has to be called often.
%They just exist to improce speed
AdjacencySmallList=[];
FacelistLocalList=[];


%A List of alle faces of all volumes 
AllFaces=zeros(MaxFacesPerVolume,MaxFaceSize,AmountOfVolumes);

%A List of alle faces of all volumes  with just local information
AllFacesSmall=zeros(MaxFacesPerVolume,MaxFaceSize,AmountOfVolumes);

%A list of all faces
FaceList=zeros(AmountOfFaces,MaxFaceSize);
counter=1;

%goes through all volumes an compute their faces
for i=1:AmountOfVolumes
    
    %the local adjacency matrix
    ALocal=Lattice(:,i);
    ALocal=reshape(ALocal,[AmountOfVertices,AmountOfVertices]);

    %The indices alias the ponts of that volume
    Indices=find(sum(ALocal));

    %The local adjacency matrix without zero lines and colmns
    ALocalSmall=ALocal(Indices,Indices);

    %Technically the same as ALocalSmall
    ALocalSmallEmbedded=zeros(MaxVerticesPerVolume);
    ALocalSmallEmbedded(1:length(ALocalSmall),1:length(ALocalSmall))=ALocalSmall;
    found=false;
    foundValue=0;

    %searches, if the structure of the adjacency matrix already was used
    [a,~,cc]=size(AdjacencySmallList);
    if a>0
    for j=1:cc
        if ~found
            if max(max(AdjacencySmallList(:,:,j)-ALocalSmallEmbedded))==0
                found=true;
                foundValue=j;
            end
        end
    end
    end

    %if so, computation is much easier
    if found

        %get the local face list
        FacelistLocal=FacelistLocalList(:,:,foundValue);
        FacelistLocal=FacelistLocal(sum(FacelistLocal,2)>0,:);
    else
        %compute the local face list
        FacelistLocal=computeFaceMatrix(ALocalSmall,false);

        %Technically the same as FacelistLocal
        FacelistLocalEmbedded=zeros(MaxFacesPerVolume,4);
        FacelistLocalEmbedded(1:length(FacelistLocal),:)=FacelistLocal;

        %if there was no entry before...
        if isempty(FacelistLocalList)

            %...create a new entry
            FacelistLocalList=FacelistLocalEmbedded;
        else

            %otherwise add the entry
            FacelistLocalList=cat(3,FacelistLocalList,FacelistLocalEmbedded);
        end

        %if there was no entry before...
        if isempty(AdjacencySmallList)

            %...create a new entry
            AdjacencySmallList=ALocalSmallEmbedded;
        else

            %otherwise add the entry
            AdjacencySmallList=cat(3,AdjacencySmallList,ALocalSmallEmbedded);
        end 
        
        

    end
    
    %The size of the local face list
    [c,d]=size(FacelistLocal);

    %Add the values 
    AllFacesSmall(1:c,1:d,i)=FacelistLocal;

    %replace the local values by the indices
    [a,b]=size(FacelistLocal);
    for j=1:a
        for k=1:b
            if FacelistLocal(j,k)>0
                FacelistLocal(j,k)=Indices(FacelistLocal(j,k));
            end
        end
    end

    %add the values to the global face list
    FaceList(counter:counter+a-1,1:b)=FacelistLocal;
    counter=counter+a;

    [c,d]=size(FacelistLocal);

    %add the values to the global face list
    AllFaces(1:c,1:d,i)=FacelistLocal;
end

%remove double faces
[FaceList,~,InnerInformation]=unique(FaceList,'rows');

[a,~]=size(FaceList);

IsInnerFace=zeros(a,1);
for i=1:a
    if nnz(InnerInformation==i)==2
        IsInnerFace(i)=1;
    end
end

IsInnerVertex=zeros(AmountOfVertices,1);
for i=1:AmountOfVertices
    [row,~]=find(FaceList==i);
    if sum(IsInnerFace(row))==length(row) && ~isempty(row)
        IsInnerVertex(i)=1;
    end
end


%-------------------------------------------------------------------------
%       Compute the new control points
%-------------------------------------------------------------------------


[S,LatticeNew]=computeTriCubicLatticeAndSubdivisionMatrix(Lattice,find(IsInnerVertex));
ControlpointsNew=S*Controlpoints;

end






%Also a subfunction in computeTriCubicSubdivisionMatrixBig 
function [S,LatticeNew] = computeTriCubicLatticeAndSubdivisionMatrix(Lattice,RefinementAreas)
%COMPUTETRICUBICLATTICEANDSUBDIVISIONMATRIX computes the subdivision matrix
%and the new lattice for a refinement in the given refinement area.
%
% Input:    Lattice:            The given lattice. Each colums is the 
%                               adjacency matrix of a cube
%           RefinementAreas:    Numbers of points arround which a
%                               refinement should be created. These points
%                               are the center of a small subdivision
%                               matrix
%
%Output:    S:                  Subdivision matrix of all points of all
%                               areas, which should refined.
%           Lattice:            The lattice of the refined area.


%The Amount of control points and the amount of cubes
[AmountOfVertices,AmountOfVolumes]=size(Lattice);
AmountOfVertices=sqrt(AmountOfVertices);


%The adjacency matrix of the whole structure
ATotal=sum(Lattice,2);
ATotal=reshape(ATotal,[AmountOfVertices,AmountOfVertices]);

%determines, which point is a corner of which cube
ActivePoints=zeros(8,AmountOfVolumes);

%fill the list of active points
for i=1:AmountOfVolumes
    LocalCube=reshape(Lattice(:,i),[AmountOfVertices,AmountOfVertices]);
    LocalCube=sum(LocalCube);
   
    Entries=find(LocalCube)';
    ActivePoints(1:length(Entries),i)=Entries;
end


%compute the refinement matrix of the regular area
F=computePrismFaceMatrix(4);
FA=FaceToAdjacencyMatrix(F);
[SS2Alt2Reg,AdjacencyMatrixCubicRegular] = computeTriCubicSubdivisionMatrix(FA,'PreventInputCheck');




%compute the regular graph
GA=graph(FA);

%compute the regular face matrix
FaceMatrixReg=computeFaceMatrix(FA,false);


%-----------------------------------------------------------------
% Precompute the necessary semi-irregular subdivision matrices
%-----------------------------------------------------------------

%Get the amount of cubes arround each vertex
[~,~,idx] = unique(reshape(ActivePoints,[8*AmountOfVolumes,1]));
AmountOfCubesArroundVertex = accumarray(idx(:),1);

% Get the vertices, which have more than 8 cubes
VerticesWithManyCubes=find(AmountOfCubesArroundVertex>8);

%The potential semi-irregular matrixes to compute
PotentialFaces=zeros(max(VerticesWithManyCubes),1);

%3,4,5 are computed as standard
PotentialFaces(3:5)=1;

%goes through all vertices with many cubes
for i = 1: length(VerticesWithManyCubes)

    %get the dual structure
    D=computeDualStructureAroundVertex(ActivePoints,VerticesWithManyCubes(i));

    %The structure of D has to be valid 
    if is3Connected(D) && PlanarityTest(D)
        %get the faces of this structure
        FZentral=computeFaceMatrix(D,false);
        
        %The amount of faces of dhe dual structore
        [AmountOfFacesLocal,~]=size(FZentral);
        
        %get the valence of the faces
        FaceNumbersLocal=zeros(AmountOfFacesLocal,1);
        for j=1:AmountOfFacesLocal
            FaceNumbersLocal(j)=nnz(FZentral(j,:));
        end
        
        %unique the valences
        FaceNumbersLocal=unique(FaceNumbersLocal);
    
        %add the number to the potential face list
        PotentialFaces(FaceNumbersLocal)=1;
    end

end

%The values of semi-irregular subdivision matrices, which has to compute
FaceNumbersZentral=find(PotentialFaces);

%The list of subdivision matrices
SSemi=cell(max(FaceNumbersZentral),1);

%The graphs of the structure of the subdivision matrices
GASemi=cell(max(FaceNumbersZentral),1);

%The face matrix of the subdivision matrices
FaceMatrixSemi=cell(max(FaceNumbersZentral),1);

%the adjacency matrix of the subdivision matrices
AdjacencyMatrixSemi=cell(max(FaceNumbersZentral),1);

%goes through all potential valences to compute 
for i=1:max(FaceNumbersZentral)

    %if the valence should be computed
    if ismember(i,FaceNumbersZentral)
    
        %compute the face matrix
        FDualLocal=computePrismFaceMatrix(i);

        %compute the adjacency matrix
        ADualLocal=FaceToAdjacencyMatrix(FDualLocal);

        %compute the subdivision matrix
        [SLocal,AdjacencyMatrixLocal] = computeTriCubicSubdivisionMatrix(ADualLocal,'PreventInputCheck');
        

        %compute the graph of the dual structure
        GALocal=graph(ADualLocal);

        %compute the face matrix of the dual structure
        FaceMatrixLocal=computeFaceMatrix(ADualLocal,false);

        %set the values
        SSemi(i)={SLocal};
        GASemi(i)={GALocal};
        FaceMatrixSemi(i)={FaceMatrixLocal};
        FaceMatrixReg=computeFaceMatrix(FA,false);
        AdjacencyMatrixSemi(i)={AdjacencyMatrixLocal};
    end
end


%-----------------------------------------------------------------
% Compute all faces of all cubes
%-----------------------------------------------------------------

%the total amount of faces
AmountOfFaces=sum((sum(Lattice)/2)/2+2-2);

%the size of the faces (as all volumes are cubes it is 4)
MaxFaceSize=4;


%local variables
AdjacencySmallList=[];
FacelistLocalList=[];


%The list of all faces
FaceList=zeros(AmountOfFaces,MaxFaceSize);
counter=1;

for i=1:AmountOfVolumes

    %the local adjacency matrix of the cube
    ALocalBig=Lattice(:,i);
    ALocalBig=reshape(ALocalBig,[AmountOfVertices,AmountOfVertices]);
    Indices=find(sum(ALocalBig));

    %a small version of the adjacency matrix without zero rows and columns
    ALocalSmall=full(ALocalBig(Indices,Indices));

    %tries to find this permutation of the cubic adjacency matrix
    found=false;
    foundValue=0;
    [a,~,cc]=size(AdjacencySmallList);
    if a>0
    for j=1:cc
        if ~found
            if max(max(AdjacencySmallList(:,:,j)-ALocalSmall))==0
                found=true;
                foundValue=j;
            end
        end
    end
    end

    %if so, faces can be set easy
    if found
        FacelistLocal=FacelistLocalList(:,:,foundValue);
        FacelistLocal=FacelistLocal(sum(FacelistLocal,2)>0,:);

    %if not, faces has to be computed
    else
        %compute the face matrix
        FacelistLocal=computeFaceMatrix(ALocalSmall,false);

        %set it to the list
        if isempty(FacelistLocalList)
            FacelistLocalList=FacelistLocal;
        else
            FacelistLocalList=cat(3,FacelistLocalList,FacelistLocal);
        end

        %set the small adjacency matrix to the list
        if isempty(AdjacencySmallList)
            AdjacencySmallList=ALocalSmall;
        else
            AdjacencySmallList=cat(3,AdjacencySmallList,ALocalSmall);
        end 
        
        

    end
    
    %Rearange the indices of the faces to the global context
    [a,b]=size(FacelistLocal);
    for j=1:a
        for k=1:b
            if FacelistLocal(j,k)>0
                FacelistLocal(j,k)=Indices(FacelistLocal(j,k));
            end
        end
    end

    %set the face list to the global list
    FaceList(counter:counter+a-1,1:b)=FacelistLocal;

    counter=counter+a;
end

%uniques the faces
FaceList=unique(FaceList,'rows');

%sort the entries for unification
FaceList=sort(FaceList,2);

%compute an unique indice for each face
FaceListNumber=(FaceList(:,1)-1)*AmountOfVertices^3+(FaceList(:,2)-1)*AmountOfVertices^2+(FaceList(:,3)-1)*AmountOfVertices+FaceList(:,4);

%the maximum amout of new points out of the face (as it is computed by the
%"edge rule" it is called edge point)
AmountOfNewEdgePoints=length(FaceListNumber);

%-----------------------------------------------------------------
% Compute additional variables
%-----------------------------------------------------------------

%the maximal amount of new volume points (is equal to the amount of points.
%As it is the "volume rule" it is called volume point)
[AmountOfNewVolumePoints,~]=size(ATotal);

%graph of the whole structure
GraphA=graph(ATotal);

%all edges of the whole structure
EdgeList=GraphA.Edges.EndNodes;

%unique number for each edge
EdgeValueNumber=(EdgeList(:,1)-1)*AmountOfVertices+EdgeList(:,2);

%the maximum amout of new points out of the edge (as it is computed by the
%"face rule" it is called face point)
AmountOfNewFacePoints=length(EdgeList);

%the maximum amout of new points out of the points (as it is computed by the
%"volume rule" it is called point point)
[~,AmountOfNewPointPoints]=size(Lattice);


%the maount of maximal new points in total
AmointOfNewPointsTotal=AmountOfNewPointPoints+AmountOfNewEdgePoints+ AmountOfNewFacePoints+AmountOfNewVolumePoints;

%determines, if for teh corresponding point there was an subdivision entry
%created
EntryCreated=zeros(AmountOfNewPointPoints+AmountOfNewEdgePoints+ AmountOfNewFacePoints+AmountOfNewVolumePoints,1);

%determines, if a refinement area already was visited
RefinementAreaVisited=zeros(1,length(RefinementAreas));

%counts the new amount of cubes 
newCubeCounter=1;


%initializing local varibales
rowAdjacency1=zeros(10000000,1);
rowAdjacency2=zeros(10000000,1);
colAdjacency=zeros(10000000,1);
entryAdjacency=zeros(10000000,1);
AdjacencyCounter=0;


%initializing local varibales
rowS=zeros(10000000,1);
colS=zeros(10000000,1);
entryS=zeros(10000000,1);
SCounter=0;


%-----------------------------------------------------------------
% the main refinement process
%-----------------------------------------------------------------


%p means priority. the first run is for regular refinement areas, the
%second one for semi-irregular and the last one for irregular. rules might
%be overrided by rules with bigger p
for p=1:3

    %goes through all refinement areas
    for i=1:length(RefinementAreas)
        
        %resets the adjacency current variable
        AdjacencyCurrent=[];

        %just goes in, if the refinement area was not visidet before (saves
        %time as things do not has to checked again for bigger p's
        if RefinementAreaVisited(i)==0
    
            %determines, if the refinement area should be refined
            refine=false;

            %the current refinement area
            j=RefinementAreas(i);
            
            %get all cubes, which have the vertex j in common
            [~,CubesHavingVertex]=find(ActivePoints==j);
        
            %dual structure arround j
            D=computeDualStructureAroundVertex(ActivePoints,j);

            %graph of the duial structure
            DGraph=graph(D);
        
            %--------------------------------------------------------------
            % checks several cases to determine, if the p is right for the
            % current structure
            %--------------------------------------------------------------

            % Area arround the vertex micht be cubic
            if length(D)==8

                %checks, if the area is cubic
                P=isomorphism(GA,DGraph);

                %If not it as an area for p==3
                if isempty(P)
                    if p==3
                        %compute subdivision matrix, adjacency matrix of
                        %the cubes and face matrix
                        [SLocal,AdjacencyCurrent]=computeTriCubicSubdivisionMatrix(D,'PreventInputCheck');
                        FaceMatrixLocal=computeFaceMatrix(D,false);
                        refine=true;
                    end
                %If so it as an area for p==1
                else
                    if p == 1
                        %get the subdivision matrix, adjacency matrix of
                        %the cubes and face matrix
                        SLocal=SS2Alt2Reg;
                        AdjacencyCurrent=AdjacencyMatrixCubicRegular;
                        CubesHavingVertex=CubesHavingVertex(P);
                        DGraph=GA;
                        FaceMatrixLocal=FaceMatrixReg;
                        refine=true;
                    end
                end
            else
                %only possible, if p~=1
                if p>1

                    %check, if a prism is possible
                    if mod(length(D),2)==0 && length(D)<=2*max(FaceNumbersZentral)

                        %checks, if the dual structure is isomorphic to a
                        %prism
                        FaceNumber=length(D)/2;
                        GLocal=GASemi{FaceNumber};
                        if FaceNumber > 2
                            P=isomorphism(GLocal,DGraph);
                        else
                            P=[];
                        end
                        
                        %if it is not isomorphic, it is an area for p==3
                        if isempty(P)
                            if p== 3
                                %compute subdivision matrix, adjacency matrix of
                                %the cubes and face matrix
                                [SLocal,AdjacencyCurrent]=computeTriCubicSubdivisionMatrix(D,'PreventInputCheck');
                                FaceMatrixLocal=computeFaceMatrix(D,false);
                                refine=true;
                            end
                        %if it is  isomorphic, it is an area for p==2
                        else
                            if p == 2
                                %get the subdivision matrix, adjacency matrix of
                                %the cubes and face matrix
                                SLocal=SSemi{FaceNumber};
                                AdjacencyCurrent=AdjacencyMatrixSemi{FaceNumber};
                                CubesHavingVertex=CubesHavingVertex(P);
                                DGraph=GLocal;
                                FaceMatrixLocal=FaceMatrixSemi{FaceNumber};
                                refine=true;
                            end
                        end
  
                    %if not it ist as an area for p==3
                    else
                        if p == 3
                            %compute subdivision matrix, adjacency matrix of
                            %the cubes and face matrix
                            [SLocal,AdjacencyCurrent]=computeTriCubicSubdivisionMatrix(D,'PreventInputCheck');
                            FaceMatrixLocal=computeFaceMatrix(D,false);
                            refine=true;
                        end
                    end
                end
            end
        
            %if the area should refined
            if refine

                %--------------------------------------------------------------
                % compute the indices for S and Lattice
                %--------------------------------------------------------------
        
                %compute the old and new point point (to identify the right 
                %indices for S and the lattice
                PointPointOld=zeros(1,length(CubesHavingVertex));
                for k=1:length(CubesHavingVertex)
                    AdjacencyCube=Lattice(:,CubesHavingVertex(k));
                    AdjacencyCube=reshape(AdjacencyCube,[AmountOfVertices,AmountOfVertices]);
                    AdjacencyCube=AdjacencyCube(ActivePoints(:,CubesHavingVertex(k)),ActivePoints(:,CubesHavingVertex(k)));
                    DistanceMatrix=AdjacencyCube+AdjacencyCube^2;
                    CentralPointInCube=(ActivePoints(:,CubesHavingVertex(k))==j);
                    PointPointInAdjacency=(DistanceMatrix(CentralPointInCube,:)==0);
                    PointPointOld(k)=ActivePoints(PointPointInAdjacency,CubesHavingVertex(k));
                end
            
                PointPointNew=CubesHavingVertex';
                
                %compute the old and new edge point (to identify the right 
                %indices for S and the lattice
                EdgeMatrixLocal=DGraph.Edges.EndNodes;
                EdgePointOld=zeros(1,length(EdgeMatrixLocal));
                EdgePointNew=zeros(1,length(EdgeMatrixLocal));
                for k=1:length(EdgeMatrixLocal)
                    
                    
                    EdgeLocal=EdgeMatrixLocal(k,:);
                    C1L=EdgeLocal(1);
                    C2L=EdgeLocal(2);
                    C1=CubesHavingVertex(C1L);
                    C2=CubesHavingVertex(C2L);
            
                    AdjacencyCube1=Lattice(:,C1);
                    AdjacencyCube1=reshape(AdjacencyCube1,[AmountOfVertices,AmountOfVertices]);
                    AdjacencyCube1=AdjacencyCube1(ActivePoints(:,C1),ActivePoints(:,C1));
                    DistanceMatrix=AdjacencyCube1+AdjacencyCube1^2;
                    CentralPointInCube=(ActivePoints(:,C1)==j);
                    EdgePointInAdjacency1=(DistanceMatrix(CentralPointInCube,:)==0);
                    NeighborC1L=AdjacencyCube1(EdgePointInAdjacency1,:);
                    NeighborC1=ActivePoints(logical(NeighborC1L),C1);
            
                    AdjacencyCube2=Lattice(:,C2);
                    AdjacencyCube2=reshape(AdjacencyCube2,[AmountOfVertices,AmountOfVertices]);
                    AdjacencyCube2=AdjacencyCube2(ActivePoints(:,C2),ActivePoints(:,C2));
                    DistanceMatrix=AdjacencyCube2+AdjacencyCube2^2;
                    CentralPointInCube=(ActivePoints(:,C2)==j);
                    EdgePointInAdjacency2=(DistanceMatrix(CentralPointInCube,:)==0);
                    NeighborC2L=AdjacencyCube2(EdgePointInAdjacency2,:);
                    NeighborC2=ActivePoints(logical(NeighborC2L),C2);
                    
                    EdgePointOld(k)=NeighborC2(ismember(NeighborC2,NeighborC1));
            
                    ActivePoints1=ActivePoints(:,C1);
                    ActivePoints2=ActivePoints(:,C2);
                    ActivePointsInnerFace=sort(ActivePoints1(ismember(ActivePoints1,ActivePoints2)));
                    ActivePointsInnerFaceNumber=(ActivePointsInnerFace(1)-1)*AmountOfVertices^3+(ActivePointsInnerFace(2)-1)*AmountOfVertices^2+(ActivePointsInnerFace(3)-1)*AmountOfVertices+ActivePointsInnerFace(4);
                    
                    EdgePointNew(k)=find(FaceListNumber==ActivePointsInnerFaceNumber);  
                end
            
            
                %compute the old and new face point (to identify the right 
                %indices for S and the lattice
                [a,~]=size(FaceMatrixLocal);
                FacePointOld=zeros(1,a);
                FacePointNew=zeros(1,a);
                for k=1:a
                    LocalFace=FaceMatrixLocal(k,:);
                    n=nnz(LocalFace);
                    LocalFace=LocalFace(1:n);
                    ActivePointsLocal=ActivePoints(:,CubesHavingVertex(LocalFace(1)));
                    for l=2:n
                        ActivePointsLocal=ActivePointsLocal(ismember(ActivePointsLocal,ActivePoints(:,CubesHavingVertex(LocalFace(l)))));
                    end
                    if ActivePointsLocal(1)==j
                        FacePointOld(k)=ActivePointsLocal(2);
                    else
                        FacePointOld(k)=ActivePointsLocal(1);
                    end
                    ActivePointsLocal=sort(ActivePointsLocal);
                    ValueNumber=(ActivePointsLocal(1)-1)*AmountOfVertices+ActivePointsLocal(2);
                    FacePointNew(k)=find(EdgeValueNumber==ValueNumber);
            
                end
            
                %compute the old and new volume point (to identify the right 
                %indices for S and the lattice
                VolumePointOld=j;
                VolumePointNew=j;

                %----------------------------------------------------------
                % set the indices and values for the new lattice
                %----------------------------------------------------------
            
                [~,~,AmountOfNewCubes]=size(AdjacencyCurrent);
        
                %goes through all cubes
                for ii=1:AmountOfNewCubes
        
                    AdjacencyLocalSmall=AdjacencyCurrent(:,:,ii);
        
                    %set the row indices
                    row1=repmat([PointPointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints+AmountOfNewEdgePoints, ...
                    EdgePointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints, ...
                    FacePointNew+AmountOfNewVolumePoints, ...
                    VolumePointNew]',length(AmointOfNewPointsTotal*([PointPointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints+AmountOfNewEdgePoints, ...
                    EdgePointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints, ...
                    FacePointNew+AmountOfNewVolumePoints, ...
                    VolumePointNew]-1)),1);
                    
                    
                    
                    row2=repelem([PointPointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints+AmountOfNewEdgePoints, ...
                    EdgePointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints, ...
                    FacePointNew+AmountOfNewVolumePoints, ...
                    VolumePointNew]',length(AmointOfNewPointsTotal*([PointPointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints+AmountOfNewEdgePoints, ...
                    EdgePointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints, ...
                    FacePointNew+AmountOfNewVolumePoints, ...
                    VolumePointNew]-1)),1);

                    %set the col indces 
                    col=repelem(newCubeCounter,length(row1),1);
                    entry=reshape(AdjacencyLocalSmall,[length(AdjacencyLocalSmall)^2,1]);
        
                    %allocate space if necessary
                    maxNewEntries=length(row1);
                    if maxNewEntries+AdjacencyCounter>length(rowAdjacency1)
                        rowAdjacency1T=rowAdjacency1;
                        rowAdjacency2T=rowAdjacency2;
                        colAdjacencyT=colAdjacency;
                        entryAdjacencyT=entryAdjacency;
                        
        
                        rowAdjacency1=zeros(length(rowAdjacency1)+10000000+maxNewEntries,1);
                        rowAdjacency2=zeros(length(rowAdjacency2)+10000000+maxNewEntries,1);
                        colAdjacency=zeros(length(colAdjacency)+10000000+maxNewEntries,1);
                        entryAdjacency=zeros(length(entryAdjacency)+10000000+maxNewEntries,1);
                        
        
                        rowAdjacency1(1:length(rowAdjacency1T))=rowAdjacency1T;
                        rowAdjacency2(1:length(rowAdjacency2T))=rowAdjacency2T;
                        colAdjacency(1:length(colAdjacencyT))=colAdjacencyT;
                        entryAdjacency(1:length(entryAdjacencyT))=entryAdjacencyT;
                        
        
                    end
        
                    %set the global values
                    rowAdjacency1(AdjacencyCounter+1:AdjacencyCounter+length(row1))=row1;
                    rowAdjacency2(AdjacencyCounter+1:AdjacencyCounter+length(row2))=row2;
                    colAdjacency(AdjacencyCounter+1:AdjacencyCounter+length(col))=col;
                    entryAdjacency(AdjacencyCounter+1:AdjacencyCounter+length(entry))=entry;
                    
                    AdjacencyCounter=AdjacencyCounter+maxNewEntries;
                    newCubeCounter=newCubeCounter+1;
                end
        
                %----------------------------------------------------------
                % set the indices and values for S
                %----------------------------------------------------------
        
                [a,b]=size(SLocal');
                maxNewEntries=a*b;

                %allocate space if necessary
                if maxNewEntries+SCounter>length(rowS)
                    rowST=rowS;
                    colST=colS;
                    entryST=entryS;
                    
        
                    rowS=zeros(length(rowS)+10000000+maxNewEntries,1);
                    colS=zeros(length(colS)+10000000+maxNewEntries,1);
                    entryS=zeros(length(entryS)+10000000+maxNewEntries,1);
                    
        
                    rowS(1:length(rowST))=rowST;
                    colS(1:length(colST))=colST;
                    entryS(1:length(entryST))=entryST;
                    
        
                end
        
                %compute the local indices
                IndicesSRow=repmat([PointPointOld,EdgePointOld,FacePointOld,VolumePointOld],1,b);
                IndicesSCol=repelem([PointPointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints+AmountOfNewEdgePoints, ...
                    EdgePointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints, ...
                    FacePointNew+AmountOfNewVolumePoints, ...
                    VolumePointNew],a);
                IndicesSEntry=reshape(SLocal',[a*b,1]);
        
        
                %set the values
                rowS(SCounter+1:SCounter+length(IndicesSRow))=IndicesSRow;
                colS(SCounter+1:SCounter+length(IndicesSCol))=IndicesSCol;
                entryS(SCounter+1:SCounter+length(IndicesSEntry))=IndicesSEntry;
                
                SCounter=SCounter+maxNewEntries;
                  
                %set the information, if an entry was created for that
                %points
                EntryCreated([PointPointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints+AmountOfNewEdgePoints, ...
                EdgePointNew+AmountOfNewVolumePoints+AmountOfNewFacePoints, ...
                FacePointNew+AmountOfNewVolumePoints, ...
                VolumePointNew])=1;
                RefinementAreaVisited(i)=1;
         
            end
        end
    end

end

%-----------------------------------------------------------------
% Set S
%-----------------------------------------------------------------

%reduce the variables (if too much space as allocated)
rowS=rowS(1:SCounter);
colS=colS(1:SCounter);
entryS=entryS(1:SCounter);

%remove empty entries
rowS(entryS==0)=[];
colS(entryS==0)=[];
entryS(entryS==0)=[];

rowS=flipud(rowS);
colS=flipud(colS);
entryS=flipud(entryS);

%put the indices in the right shape (as several rows are empty)
UniqueValues=rowS+colS*(length(EntryCreated)-1);
[~,bb,~]=unique(UniqueValues);
rowS=rowS(bb);
colS=colS(bb);
entryS=entryS(bb);

%create S
S=sparse(rowS,colS,entryS,AmountOfVertices,length(EntryCreated));
S=S(:,logical(EntryCreated));
S=S';


%-----------------------------------------------------------------
% Set Lattice
%-----------------------------------------------------------------

%reduce the variables (if too much space as allocated)
rowAdjacency1=rowAdjacency1(1:AdjacencyCounter);
rowAdjacency2=rowAdjacency2(1:AdjacencyCounter);
colAdjacency=colAdjacency(1:AdjacencyCounter);
entryAdjacency=entryAdjacency(1:AdjacencyCounter);


[AmountOfNewPoints,~]=size(S);

%get the lines of entries, which are really created/refined
row=find(EntryCreated);

[~,rowAdjacency1]=ismember(rowAdjacency1,row);
[~,rowAdjacency2]=ismember(rowAdjacency2,row);

%compute the row indices
rowAdjacency=rowAdjacency1+length(row)*(rowAdjacency2-1);

%remove empty entries
rowAdjacency(entryAdjacency==0)=[];
colAdjacency(entryAdjacency==0)=[];
entryAdjacency(entryAdjacency==0)=[];

%create the new lattice
LatticeNew=sparse(rowAdjacency,colAdjacency,entryAdjacency,AmountOfNewPoints^2,max(colAdjacency));

end





%Also a subfunction in computeTriCubicSubdivisionMatrixBig 
function [DualAdjacency] = computeDualStructureAroundVertex(ActivePoints,currentVertex)
% COMPUTEDUALSTRUCTUREAROUNDVERTEX computes the structure of volumes around
% a vertex in the cubic case. It computes the polyhedral structure the
% cubes have.
%
% Input:    ActivePoints:   The points of the volumes. Each colums stands 
%                           for one cube
%           currentVertex:  The vertex around which the structure should be
%                           computed
%
% Output:   Dual Adjacency: The adjacency matrix which codes the structure
%
% computeDualStructureAroundVertex(ActivePoints,i)  computes the structure 
% of volumes around a vertex in the cubic case. It computes the polyhedral 
% struchture the cubes have.

% searches fur cubes, which have the vertex i in his vertex set.
[~,CubesHavingVertex]=find(ActivePoints==currentVertex);

%the amount of volumes alias the amount of dual points
AmountOfDualPoints=length(CubesHavingVertex);

%initializing the dual adjacency matrix
DualAdjacency=zeros(AmountOfDualPoints);

%goes through all adjacent cubes
for i=1:AmountOfDualPoints

    %the actice points arround the cube i
    ActivePointsI=ActivePoints(ActivePoints(:,CubesHavingVertex(i))>0,CubesHavingVertex(i));
    
    %goes through all adjacent cubes
    for j=1:AmountOfDualPoints
        
        %no entries on the diagonal
        if i~=j

            %the actice points arround the cube j
            ActivePointsJ=ActivePoints(ActivePoints(:,CubesHavingVertex(j))>0,CubesHavingVertex(j));

            %if both volumes have more than 2 points in common they share a
            %face. So the dual points are neighbored
            if sum(ismember(ActivePointsI,ActivePointsJ))>2
                DualAdjacency(i,j)=1;
            end
        end
    end
end

end

