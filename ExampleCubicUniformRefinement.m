%------------------------------------------------------------------------
% Createan initial structure
%------------------------------------------------------------------------


A=FaceToAdjacencyMatrix(F);

L= refineTriCubicInitialStructure(A,F);



LNew=zeros(28,28);

LNew(2,12)=1;
LNew(12,2)=1;

LNew(2,13)=1;
LNew(13,2)=1;

LNew(2,28)=1;
LNew(28,2)=1;

LNew(7,13)=1;
LNew(13,7)=1;

LNew(7,15)=1;
LNew(15,7)=1;

LNew(7,28)=1;
LNew(28,7)=1;

LNew(3,12)=1;
LNew(12,3)=1;

LNew(3,15)=1;
LNew(15,3)=1;

LNew(3,28)=1;
LNew(28,3)=1;

LNew(23,13)=1;
LNew(13,23)=1;

LNew(23,12)=1;
LNew(12,23)=1;

LNew(23,15)=1;
LNew(15,23)=1;

LFinal=zeros(28,28,9);

LFinal(1:27,1:27,1:8)=L;

LFinal(:,:,9)=LNew;


[a,b,c]=size(LFinal);



Lattice=reshape(LFinal,[a*b,c]);






%------------------------------------------------------------------------
% Creates a List of all faces 
%------------------------------------------------------------------------




%The amount of vertices and the amount of cubes
[AmountOfVertices,AmountOfVolumes]=size(Lattice);
AmountOfVertices=sqrt(AmountOfVertices);

% The adjacency matrix of the whole structure
AdjacencyTotal=sum(Lattice,2);
AdjacencyTotal=reshape(AdjacencyTotal,[AmountOfVertices,AmountOfVertices]);

G=graph(AdjacencyTotal);

%A List of all edges of the structure
EdgeList=G.Edges.EndNodes;

%an individual number of each edge. Easier for identification
EdgeValueNumber=(EdgeList(:,1)-1)*AmountOfVertices+EdgeList(:,2);

%The amount of edges
AmountOfEdges=length(EdgeList);

%The amount of faces
AmountOfFaces=sum((sum(Lattice)/2)/2+2-2);

%As all volumes are cubes, all faces are quads
MaxFaceSize=4;

%The new amount of volumes (complicated formula for 8*Amount of volumes)
AmountOfNewVolumes=sum(sum(Lattice)/2/2+2);

%The maximum amount of faces per volume (complicated formula for 6)
MaxFacesPerVolume=max((sum(Lattice)/2)/2+2-2);

%The maximum amount of faces per volume (complicated formula for 8)
MaxVerticesPerVolume= max(sum(Lattice)/2/2+2);


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
FaceList=unique(FaceList,'rows');


%------------------------------------------------------------------------
% Creates the initial Control points
%------------------------------------------------------------------------

S=computeTriQuadraticSubdivisionMatrixV1(F,1/4);



B=S-1/2*eye(length(S));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
D=D(a)+1/2;


V=orth(real(V));


[row,col]=find(A);
EdgeList=[row,col];
EdgeList(EdgeList(:,1)>EdgeList(:,2),:)=[];
EdgeList=sortrows(EdgeList);

Controlpoints=[V;...
    V(EdgeList(:,1),:)/2+V(EdgeList(:,2),:)/2;...
    V(1,:)/5+V(2,:)/5+V(3,:)/5+V(4,:)/5+V(5,:)/5;...
    V(1,:)/5+V(2,:)/5+V(7,:)/5+V(8,:)/5+V(6,:)/5;...
    V(2,:)/3+V(3,:)/3+V(4,:)/3;...
    V(3,:)/4+V(4,:)/4+V(8,:)/4+V(7,:)/4;...
    V(4,:)/4+V(5,:)/4+V(6,:)/4+V(7,:)/4;...
    V(5,:)/3+V(1,:)/3+V(6,:)/3;...
    0,0,0;...
    2*V(2,:)/3+2*V(3,:)/3+2*V(4,:)/3];



LatticeGlobal=sum(Lattice,2);


[a,~]=size(Lattice);
LatticeGlobal=reshape(LatticeGlobal,[sqrt(a),sqrt(a)]);


G=graph(LatticeGlobal);
Edges=G.Edges.EndNodes;

ControlpointsVertices=Controlpoints;
ControlpointsEdges=1/2*Controlpoints(Edges(:,1),:)+1/2*Controlpoints(Edges(:,2),:);
ControlpointsFaces=1/4*Controlpoints(FaceList(:,1),:)+1/4*Controlpoints(FaceList(:,2),:)+1/4*Controlpoints(FaceList(:,3),:)+1/4*Controlpoints(FaceList(:,4),:);

[a,b]=size(Lattice);
IndicesVolumes=zeros(b,8);
for i=1:b
    LatticeLocal=Lattice(:,i);
    LatticeLocal=reshape(LatticeLocal,[sqrt(a),sqrt(a)]);
    [row,~]=find(LatticeLocal);
    row=unique(row);
    IndicesVolumes(i,:)=row';
end

ControlpointsVolumes=1/8*Controlpoints(IndicesVolumes(:,1),:)+1/8*Controlpoints(IndicesVolumes(:,2),:)+1/8*Controlpoints(IndicesVolumes(:,3),:)+1/8*Controlpoints(IndicesVolumes(:,4),:)+1/8*Controlpoints(IndicesVolumes(:,5),:)+1/8*Controlpoints(IndicesVolumes(:,6),:)+1/8*Controlpoints(IndicesVolumes(:,7),:)+1/8*Controlpoints(IndicesVolumes(:,8),:);

Controlpoints=[ControlpointsVertices;ControlpointsEdges;ControlpointsFaces;ControlpointsVolumes];

%------------------------------------------------------------------------
% Creates the initial lattice
%------------------------------------------------------------------------


Lattice=RefineTriCubicLatticeStructure(Lattice);


%------------------------------------------------------------------------
% Refines uniform until Spline evaluation is possible
%------------------------------------------------------------------------

[Controlpoints,Lattice] = RefineTriCubicLatticeUniform(Controlpoints,Lattice);
[Controlpoints,Lattice] = RefineTriCubicLatticeUniform(Controlpoints,Lattice);

%------------------------------------------------------------------------
% creates the figures
%------------------------------------------------------------------------

figure
plotTriCubicBSplineLattice(Controlpoints,Lattice)
title('start subdivision volume')
view(-108,27);

[Controlpoints,Lattice] = RefineTriCubicLatticeUniform(Controlpoints,Lattice);
figure
plotTriCubicBSplineLattice(Controlpoints,Lattice)
title('subdivision volume after refinement')
view(-108,27);

%-------------------------------------------------------------------------
%       Local Help Functions
%-------------------------------------------------------------------------
function [LatticeNew,ActivePoints] = RefineTriCubicLatticeStructure(Lattice)
% REFINETRUCUBICLATTICESTRUCTURE refines the given lattice globally. Creates
% for each cube 8 new cubes.
%
% Input:    Lattice:        The lattice. Each column stands for one cube.
% Output:   LatticeNew:     The new lattice
%           ActivePoints:   The point indices on each cube.
%
% RefineTriCubicLatticeStructure(Lattice) refines the given lattice 
% globally. Creates for each cube 8 new cubes.

%The amount of vertices and the amount of cubes
[AmountOfVertices,AmountOfVolumes]=size(Lattice);
AmountOfVertices=sqrt(AmountOfVertices);

% The adjacency matrix of the whole structure
AdjacencyTotal=sum(Lattice,2);
AdjacencyTotal=reshape(AdjacencyTotal,[AmountOfVertices,AmountOfVertices]);

G=graph(AdjacencyTotal);

%A List of all edges of the structure
EdgeList=G.Edges.EndNodes;

%an individual number of each edge. Easier for identification
EdgeValueNumber=(EdgeList(:,1)-1)*AmountOfVertices+EdgeList(:,2);

%The amount of edges
AmountOfEdges=length(EdgeList);

%The amount of faces
AmountOfFaces=sum((sum(Lattice)/2)/2+2-2);

%As all volumes are cubes, all faces are quads
MaxFaceSize=4;

%The new amount of volumes (complicated formula for 8*Amount of volumes)
AmountOfNewVolumes=sum(sum(Lattice)/2/2+2);

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
FaceList=unique(FaceList,'rows');

%give each face a unique number (improving speed)
FaceListNumber=(FaceList(:,1)-1)*AmountOfVertices^3+(FaceList(:,2)-1)*AmountOfVertices^2+(FaceList(:,3)-1)*AmountOfVertices+FaceList(:,4);

%the amount of all faces
[AmountOfFaces,~]=size(FaceList);

%the amount of new points
TotalNewPoints=AmountOfVertices+AmountOfEdges+AmountOfFaces+AmountOfVolumes;

%-------------------------------------------------------------------------
%       Compute the adjacency matrix and the active points
%-------------------------------------------------------------------------

counter=0;
IndicesRowTotal=zeros(1,24*AmountOfNewVolumes);
IndicesColTotal=zeros(1,24*AmountOfNewVolumes);
ActivePoints=zeros(8,AmountOfNewVolumes);
ActivePointsCounter=0;
IndicesCounter=0;

%goes through each volume and divide it into 8 cubes
for i=1:AmountOfVolumes

    
    %compute local information

    %local adjacency matrix
    ALocal=Lattice(:,i);
    ALocal=reshape(ALocal,[AmountOfVertices,AmountOfVertices]);

    Indices=find(sum(ALocal));

    %small version
    ALocalSmall=ALocal(Indices,Indices);

    %local face version
    AllFacesSmallLocal=AllFacesSmall(:,:,i);
    AllFacesSmallLocal=AllFacesSmallLocal(sum(AllFacesSmallLocal,2)>0,:);

    %adjacency matrix of refinement
    ANewLocalSmall = refineTriCubicInitialStructure(ALocalSmall,AllFacesSmallLocal);

    %the amount of new points and volumes
    [NewPoints,~,NewVolumes]=size(ANewLocalSmall);

  
    %The local face information
    FLocal=AllFaces(:,:,i);
    FLocal=FLocal(sum(FLocal,2)>0,:);
    
    %Local edges
    [row,col]=find(ALocalSmall);
    EdgeListLocal=[row,col];
    EdgeListLocal(EdgeListLocal(:,1)>EdgeListLocal(:,2),:)=[];

    EdgeListLocal=sortrows(EdgeListLocal);

    
    %information of point indices
    IndicesToBigMatrixV=Indices;

    %global new edges
    EdgeListLocal=Indices(EdgeListLocal);

    %edge numbers
    EdgeNumberLocal=(EdgeListLocal(:,1)-1)*AmountOfVertices+EdgeListLocal(:,2);

    %the concrete edge numbers for which new points are created
    [~,Numbers]=ismember(EdgeNumberLocal,EdgeValueNumber);

    %information of point and edge indices
    IndicesToBigMatrixVE=[IndicesToBigMatrixV,AmountOfVertices+Numbers'];

   
    %get the numbers of the local faces for which new points are created
    FaceListNumberLocal=(FLocal(:,1)-1)*AmountOfVertices^3+(FLocal(:,2)-1)*AmountOfVertices^2+(FLocal(:,3)-1)*AmountOfVertices+FLocal(:,4);

    [~,Numbers]=ismember(FaceListNumberLocal,FaceListNumber);

    %information of point and edge and face indices
    IndicesToBigMatrixVEF=[IndicesToBigMatrixVE,AmountOfVertices+AmountOfEdges+Numbers'];

    %information of point and edge and face and volume indices
    IndicesToBigMatrixVEFP=[IndicesToBigMatrixVEF,AmountOfVertices+AmountOfEdges+AmountOfFaces+i];


    %The concrete indices for the new adjacency matrix
    IndicesToBigMatrixCalculated=((IndicesToBigMatrixVEFP'-1)*TotalNewPoints)+IndicesToBigMatrixVEFP;
    IndicesToBigMatrixCalculated=reshape(IndicesToBigMatrixCalculated,[NewPoints^2,1]);

    IndicesToBigMatrixCalculatedCol=ones(length(IndicesToBigMatrixCalculated),1)*(1:NewVolumes);
    IndicesToBigMatrixCalculatedCol=reshape(IndicesToBigMatrixCalculatedCol,[length(IndicesToBigMatrixCalculated)*NewVolumes,1]);

    IndicesToBigMatrixCalculated=repmat(IndicesToBigMatrixCalculated,NewVolumes,1);
    
    %the new adjacency matrix
    ANewLocal=sparse(IndicesToBigMatrixCalculated,IndicesToBigMatrixCalculatedCol,reshape(ANewLocalSmall,[NewPoints^2*NewVolumes,1]),TotalNewPoints^2,NewVolumes);

    %coded information about the new poitns on the new volumes
    ANewLocal2=sparse(mod(IndicesToBigMatrixCalculated-1,TotalNewPoints)+1,IndicesToBigMatrixCalculatedCol,reshape(ANewLocalSmall,[NewPoints^2*NewVolumes,1]),TotalNewPoints,NewVolumes);

    %translate the new adjacency matrix to row and col information
    [row,col]=find(ANewLocal);

    %set the indices for the final matrix
    IndicesRowTotal(:,IndicesCounter+1:IndicesCounter+24*8)=row;
    IndicesColTotal(:,IndicesCounter+1:IndicesCounter+24*8)=col+counter;
    IndicesCounter=IndicesCounter+24*8;
    
    counter=counter+NewVolumes;

    %set the active points
    [row,~]=find(ANewLocal2);
    ActivePoints(:,ActivePointsCounter+1:ActivePointsCounter+8)=reshape(row,[8,8]);
    ActivePointsCounter=ActivePointsCounter+8;
   
end

%set the final matrix
LatticeNew=sparse(IndicesRowTotal,IndicesColTotal,ones(length(IndicesRowTotal),1),TotalNewPoints^2,AmountOfNewVolumes);



end





function [Lattice] = refineTriCubicInitialStructure(AdjacencyMatrix,FaceMatrix)
% REFINETRICUBICINITIALSTRUCTURE refines the structure of an initial 
% tri-cubic element to an adjacency matrix structure, which codes each cube
%
% Input:    AdjacencyMatrix, FaceMatrix
% Output:   LatticeStructure of cubes

%The amount of vertices
AmountOfVertices=length(AdjacencyMatrix);

%the amount of edges
AmountOfEdges=sum(sum(AdjacencyMatrix))/2;

%if there is no face matrix given
if isempty(FaceMatrix)

    %compute the face matrix
    FaceMatrix=computeFaceMatrix(AdjacencyMatrix,'false');
end

%the amount of faces
[AmountOfFaces,~]=size(FaceMatrix);

%the amount of volumes
AmountOfVolumes=1;

%The list of edges
[row,col]=find(AdjacencyMatrix);
EdgeList=[row,col];
EdgeList(EdgeList(:,1)>EdgeList(:,2),:)=[];
EdgeList=sortrows(EdgeList);


%The amount of points the new structure has (for ever point, edge, face and
%vertex there is a new vertex)
TotalNewPoints=AmountOfVertices+AmountOfEdges+AmountOfFaces+AmountOfVolumes;

%the new structure matrix. Each (:,:,i) element is the adjacency matrix of
%one cube
Lattice=zeros(TotalNewPoints,TotalNewPoints,AmountOfVertices);

for i=1:AmountOfVertices
    [ConnectedEdges,~]=find(EdgeList==i);    

    %computes the adjacency entries for the edge of a primal point to its
    %edge points
    for j=1:length(ConnectedEdges)
        Lattice(i,AmountOfVertices+ConnectedEdges(j),i)=1;
        Lattice(AmountOfVertices+ConnectedEdges(j),i,i)=1;
    end

    

    for j=1:length(ConnectedEdges)

        %searches for the faces on that edge
        PotentialFace=zeros(2,1);
        PotentialFaceCounter=1;
        EdgeTemp=EdgeList(ConnectedEdges(j),:);
        for k=1:AmountOfFaces
            if sum(FaceMatrix(k,:)==EdgeTemp(1))+sum(FaceMatrix(k,:)==EdgeTemp(2))==2
                PotentialFace(PotentialFaceCounter)=k;
                PotentialFaceCounter=PotentialFaceCounter+1;
            end
        end


        if length(PotentialFace)==2

            %computes the adjacency entries for the edge of a primal point to its
            %face points
            Lattice(AmountOfVertices+ConnectedEdges(j),AmountOfVertices+AmountOfEdges+PotentialFace(1),i)=1;
            Lattice(AmountOfVertices+AmountOfEdges+PotentialFace(1),AmountOfVertices+ConnectedEdges(j),i)=1;

            Lattice(AmountOfVertices+ConnectedEdges(j),AmountOfVertices+AmountOfEdges+PotentialFace(2),i)=1;
            Lattice(AmountOfVertices+AmountOfEdges+PotentialFace(2),AmountOfVertices+ConnectedEdges(j),i)=1;

            %computes the adjacency entries for the face of a primal point to its
            %volume points
            Lattice(end,AmountOfVertices+AmountOfEdges+PotentialFace(1),i)=1;
            Lattice(AmountOfVertices+AmountOfEdges+PotentialFace(1),end,i)=1;

            Lattice(end,AmountOfVertices+AmountOfEdges+PotentialFace(2),i)=1;
            Lattice(AmountOfVertices+AmountOfEdges+PotentialFace(2),end,i)=1;

        else
            error('Structure was not valid.')
        end
    end

end



end