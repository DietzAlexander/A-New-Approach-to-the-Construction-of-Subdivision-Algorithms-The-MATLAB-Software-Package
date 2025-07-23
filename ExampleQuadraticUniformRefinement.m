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




[Controlpoints,Lattice] = RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,1/4);


[Controlpoints,Lattice] = RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,1/4);
figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice)
title('start subdivision volume')
view(-170,53)
[Controlpoints,Lattice] = RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,1/4);
figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice)
title('subdivision volume after refinement')
view(-170,53)

%-------------------------------------------------------------------------
%       Local Help Functions
%-------------------------------------------------------------------------


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