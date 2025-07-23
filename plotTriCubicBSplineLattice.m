function [] = plotTriCubicBSplineLattice(Controlpoints,Lattice,varargin)
% PLOTTRICUBICBSPLINELATTICE plots the plotable area of a given
% cubic B-Spline lattice.
%
% Input:    Controlpoints:  The controlpoints of the lattice
%           Lattice:        The information of geometry. Each column is an
%                           adjacency matrix of one volume of the lattice.
%
% plotTriCubicBSplineLattice(Controlpoints,Lattice) plots the plotable 
% area of a given cubic B-Spline lattice.
%
% plotTriCubicBSplineLattice( ,'PointsPerCube',n) plots the B-Spline
% lattice with nxnxn points per cube.
%
%plotTriCubicBSplineLattice( ,'ColorTable',CT) plots the B-Spline
% lattice with the colors in the CT. Each row is one color. The color is
% changed for each cube. If there are more Cubes then colors the color
% starts from the beginning.

%Input check
if max(max(abs(imag(Controlpoints)))) >0
 warning ('Controlpoints have an imaginary part. This imaginary part is ignored.')
 Controlpoints=real(Controlpoints);
end


%The amount of additional input parameters  
NumberOfAdditionalInput = nargin - 3;
PointsPerCube=7;

%Different green colors. They are used for the different cubes, which are
%plotted
ColorTable=[0.8,1,0.8;...  %1
            0.7,0.9,0.7;...%2
            0.6,0.8,0.6;...%3
            0.5,0.7,0.5;...%4
            0.4,0.6,0.4;...%5
            0.3,0.5,0.3;...%6
            0.2,0.4,0.2;...&7
            0.1,0.3,0.1;...%8
            0.0,0.2,0.0;...%9
            0,1/2,0;...    %10
            0.5,1,0.5;...  %11
            0.4,0.9,0.4;...%12
            0.3,0.8,0.3;...%13
            0.2,0.7,0.2;...%14
            0.1,0.6,0.1;...%15
            ];

i=1;
while i <= NumberOfAdditionalInput

    %Specifies whether a status of the algorithm progress should be shown
    if ischar(varargin{i}) && strcmp(varargin{i},'PointsPerCube')
        PointsPerCube=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;

    %Specifies if the central subdivision algorithm should do an input
    %check
    elseif ischar(varargin{i}) && strcmp(varargin{i},'ColorTable')
        ColorTable=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;
    else
        i=i+1;
    end
end

%----------------------------------------------------------------------
%  Compute Data Structure
%----------------------------------------------------------------------

%Amount of several objects
[AmountOfVertices,AmountOfVolumes]=size(Lattice);
AmountOfVertices=sqrt(AmountOfVertices);
AmountOfFaces=sum((sum(Lattice)/2)/2+2-2);
MaxFacesPerVolume=max((sum(Lattice)/2)/2+2-2);
MaxVerticesPerVolume= max(sum(Lattice)/2/2+2);

%Variables for the next loop
AdjacencySmallList=[];
FacelistLocalList=[];
MaxFaceSize=4;
AllFaces=zeros(MaxFacesPerVolume,MaxFaceSize,AmountOfVolumes);



FaceList=zeros(AmountOfFaces,MaxFaceSize);
counter=1;

%Goes through all cubes and computes the faces
for i=1:AmountOfVolumes

    %The local adjacency matrx
    ALocal=Lattice(:,i);
    ALocal=reshape(ALocal,[AmountOfVertices,AmountOfVertices]);

    Indices=find(sum(ALocal));
    ALocalSmall=ALocal(Indices,Indices);

    ALocalSmallEmbedded=zeros(MaxVerticesPerVolume);
    ALocalSmallEmbedded(1:length(ALocalSmall),1:length(ALocalSmall))=ALocalSmall;
    found=false;
    foundValue=0;
    [a,~,cc]=size(AdjacencySmallList);

    %If the local adjacency matrix was already there in the loop, things
    %can speed up. So check, if the local adjacency matrix was seen before
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

    %if found, get the values
    if found
        FacelistLocal=FacelistLocalList(:,:,foundValue);
        FacelistLocal=FacelistLocal(sum(FacelistLocal,2)>0,:);

    %if not, compute the values
    else

        %compute the local face list, saves it and saves the corresponding
        %adjacency matrix
        FacelistLocal=computeFaceMatrix(ALocalSmall,false);
        FacelistLocalEmbedded=zeros(MaxFacesPerVolume,4);
        FacelistLocalEmbedded(1:length(FacelistLocal),:)=FacelistLocal;
        if isempty(FacelistLocalList)
            FacelistLocalList=FacelistLocalEmbedded;
        else
            FacelistLocalList=cat(3,FacelistLocalList,FacelistLocalEmbedded);
        end

        if isempty(AdjacencySmallList)
            AdjacencySmallList=ALocalSmallEmbedded;
        else
            AdjacencySmallList=cat(3,AdjacencySmallList,ALocalSmallEmbedded);
        end 
        
        
    end
    
    %get the global values
    [a,b]=size(FacelistLocal);
    for j=1:a
        for k=1:b
            if FacelistLocal(j,k)>0
                FacelistLocal(j,k)=Indices(FacelistLocal(j,k));
            end
        end
    end

    %set the values of the faces to the list of all faces
    FaceList(counter:counter+a-1,1:b)=FacelistLocal;
    [c,d]=size(FacelistLocal);
    AllFaces(1:c,1:d,i)=FacelistLocal;

    counter=counter+a;
end

%get a unique number to each face to speed up 
FaceListNumber=(FaceList(:,1)-1)*AmountOfVertices^3+(FaceList(:,2)-1)*AmountOfVertices^2+(FaceList(:,3)-1)*AmountOfVertices+FaceList(:,4);

%compute the information, which face is an inner face
[~,~,idx] = unique(FaceListNumber);
InnerFaces = accumarray(idx(:),1);
InnerFaces=InnerFaces>1;


[FaceList,~,~]=unique(FaceList,'rows');


%----------------------------------------------------------------------
%  plot the B-Spline volumes
%----------------------------------------------------------------------

%The amount of colors
[AmountOfColors,~]=size(ColorTable);

%The evaluation points in the domain
EvaluationPoints=0:1/(PointsPerCube-1):1;

                         
%Matrix that converges polynoms to B-Splines
GGGG=[ 1/6, -1/2,  1/2, -1/6;... 
       2/3,  0,   -1,    1/2;...
       1/6,  1/2,  1/2, -1/2;...
       0,    0,    0,    1/6];



%B-Spline Evaluation points
EvaluationPoints=(EvaluationPoints').^(0:3)*GGGG';


%Evaluation matrix
M=computeEvaluationpoints3D(eye(64),EvaluationPoints,EvaluationPoints,EvaluationPoints,4);

%computes the structure of  4x4x4 controlpoints
AdjacencyEvaluation=zeros(64);
for i=1:4
    for j=1:4
        for k=1:3
            P1=(i-1)*16+(j-1)*4+k;
            P2=(i-1)*16+(j-1)*4+k+1;

            AdjacencyEvaluation(P1,P2)=1;
            AdjacencyEvaluation(P2,P1)=1;

            P1=(i-1)*16+(k-1)*4+j;
            P2=(i-1)*16+(k-1)*4+j+4;

            AdjacencyEvaluation(P1,P2)=1;
            AdjacencyEvaluation(P2,P1)=1;

            P1=(k-1)*16+(i-1)*4+j;
            P2=(k-1)*16+(i-1)*4+j+16;

            AdjacencyEvaluation(P1,P2)=1;
            AdjacencyEvaluation(P2,P1)=1;
        end
    end
end

%the corresponding graph
GEvaluation=graph(AdjacencyEvaluation);

%computes some basic variables
PointsPerCubeLattice=sum(Lattice)/2/2+2;

MaxPoints=max(PointsPerCubeLattice);

ActivePoints=zeros(MaxPoints,AmountOfVolumes);

for i=1:AmountOfVolumes
    LocalCube=reshape(Lattice(:,i),[AmountOfVertices,AmountOfVertices]);
    LocalCube=sum(LocalCube);
    ActivePoints(:,i)=find(LocalCube)';
end

[AmountOfPoints,~]=size(Controlpoints);

ColumnIndices=repmat(1:AmountOfVolumes,MaxPoints,1);

%Entry (i,j) is one, if volume(i) has point(j)
PointsInCubes=sparse(reshape(ColumnIndices,[MaxPoints*AmountOfVolumes,1]),reshape(ActivePoints,[MaxPoints*AmountOfVolumes,1]),ones(MaxPoints*AmountOfVolumes,1),AmountOfVolumes,AmountOfPoints);

[AmountOfFaces,~]=size(FaceList);

ColumnIndices=repmat((1:AmountOfFaces),4,1);

%Entry (i,j) is one, if face(i) has point(j)
PointsInFaces=sparse(reshape(ColumnIndices,[4*AmountOfFaces,1]),reshape(FaceList',[4*AmountOfFaces,1]),ones(4*AmountOfFaces,1),AmountOfFaces,AmountOfPoints);

[a,b]=size(ActivePoints);
PointsInVolume=full(sparse(reshape(ActivePoints,[a*b,1]), ones(a*b,1), ones(a*b,1)));


colorCounter=1;

%Variables for the final surf
PlotPointsAGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube);
PlotPointsBGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube);
PlotPointsCGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube);

ColorGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube,3);

NaNCounter=1;



%goes through each cube and checks, if its an inner cube. If so, the
%structure arround it can be plotted
for i=1:AmountOfVolumes
    
    
    %the points of the current cube
    LocalCubePoints=ActivePoints(:,i);

    %the cubes which have at least one point in common with the current
    %cube
    CubeGroup=find(sum(PointsInCubes(:,LocalCubePoints),2));
   
    %if there are more or less then 27, the area could not be plotted
    if length(CubeGroup)==27
        CubeGroupPoints=unique(ActivePoints(:,CubeGroup));
        valid=true;
        for j=1:length(CubeGroupPoints)
            row=PointsInFaces(:,CubeGroupPoints(j));
    
            %if all faces of the point are inner faces
            if min(InnerFaces(logical(row)))==1
                if valid
                    %checks if the inner structure consits out of 8 cubes
                    if PointsInVolume(CubeGroupPoints(j))~=8%  length(SurroundingOfJ) ~=8
                        %if not, the area could ot be plotted
                        valid=false;
                    
                    end
    
                end
            end
        end


        %if the area can be plotted
        if valid
            %gets the local adjacency matrx
            AdjacencyLocal=sum(Lattice(:,CubeGroup),2);
            AdjacencyLocal=reshape(AdjacencyLocal,[AmountOfVertices,AmountOfVertices]);

            %get the indces of the control points
            Indices=find(sum(AdjacencyLocal));
            AdjacencyLocalSmall=AdjacencyLocal(Indices,Indices);

            %get the graph structure
            GA=graph(AdjacencyLocalSmall);

            %map the structure to the upper computed one
            I=isomorphism(GEvaluation,GA);
            EvaluationNumbers=Indices(I);

            %compute the concrete plot points
            PlotPoints=M*Controlpoints(EvaluationNumbers,:);

            %devide it into x,y and z coordinates
            PlotPointsA=reshape(PlotPoints(:,1),[PointsPerCube,PointsPerCube,PointsPerCube]);
            PlotPointsB=reshape(PlotPoints(:,2),[PointsPerCube,PointsPerCube,PointsPerCube]);
            PlotPointsC=reshape(PlotPoints(:,3),[PointsPerCube,PointsPerCube,PointsPerCube]);


            %Face 1
            PlotPointsAGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsA(1,:,:),[PointsPerCube,PointsPerCube]);
            PlotPointsBGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsB(1,:,:),[PointsPerCube,PointsPerCube]);
            PlotPointsCGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsC(1,:,:),[PointsPerCube,PointsPerCube]);

            ColorGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:,:)=repmat(reshape(ColorTable(mod(colorCounter,AmountOfColors)+1,:),[1,1,3]),PointsPerCube,PointsPerCube);

            NaNCounter=NaNCounter+PointsPerCube+1;


            %Face 2
            PlotPointsAGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsA(end,:,:),[PointsPerCube,PointsPerCube]);
            PlotPointsBGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsB(end,:,:),[PointsPerCube,PointsPerCube]);
            PlotPointsCGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsC(end,:,:),[PointsPerCube,PointsPerCube]);

            ColorGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:,:)=repmat(reshape(ColorTable(mod(colorCounter,AmountOfColors)+1,:),[1,1,3]),PointsPerCube,PointsPerCube);

            NaNCounter=NaNCounter+PointsPerCube+1;



            %Face 3
            PlotPointsAGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsA(:,1,:),[PointsPerCube,PointsPerCube]);
            PlotPointsBGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsB(:,1,:),[PointsPerCube,PointsPerCube]);
            PlotPointsCGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsC(:,1,:),[PointsPerCube,PointsPerCube]);

            ColorGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:,:)=repmat(reshape(ColorTable(mod(colorCounter,AmountOfColors)+1,:),[1,1,3]),PointsPerCube,PointsPerCube);

            NaNCounter=NaNCounter+PointsPerCube+1;


            %Face 4
            PlotPointsAGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsA(:,end,:),[PointsPerCube,PointsPerCube]);
            PlotPointsBGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsB(:,end,:),[PointsPerCube,PointsPerCube]);
            PlotPointsCGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsC(:,end,:),[PointsPerCube,PointsPerCube]);

            ColorGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:,:)=repmat(reshape(ColorTable(mod(colorCounter,AmountOfColors)+1,:),[1,1,3]),PointsPerCube,PointsPerCube);

            NaNCounter=NaNCounter+PointsPerCube+1;


            %Face 5
            PlotPointsAGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsA(:,:,1),[PointsPerCube,PointsPerCube]);
            PlotPointsBGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsB(:,:,1),[PointsPerCube,PointsPerCube]);
            PlotPointsCGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsC(:,:,1),[PointsPerCube,PointsPerCube]);

            ColorGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:,:)=repmat(reshape(ColorTable(mod(colorCounter,AmountOfColors)+1,:),[1,1,3]),PointsPerCube,PointsPerCube);

            NaNCounter=NaNCounter+PointsPerCube+1;


            %Face 6
            PlotPointsAGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsA(:,:,end),[PointsPerCube,PointsPerCube]);
            PlotPointsBGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsB(:,:,end),[PointsPerCube,PointsPerCube]);
            PlotPointsCGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:)=reshape(PlotPointsC(:,:,end),[PointsPerCube,PointsPerCube]);

            ColorGlobal(NaNCounter:NaNCounter+PointsPerCube-1,:,:)=repmat(reshape(ColorTable(mod(colorCounter,AmountOfColors)+1,:),[1,1,3]),PointsPerCube,PointsPerCube);

            NaNCounter=NaNCounter+PointsPerCube+1;

            colorCounter=colorCounter+1;
            
        end
    end
end

%Plots finally the structure
surf(PlotPointsAGlobal,PlotPointsBGlobal,PlotPointsCGlobal,ColorGlobal)
axis equal
axis off


end













%Also a local function in plotTriQuadraticBSplineLattice
function [M] = computeEvaluationpoints3D(P,x,y,z,dim)
% COMPUTEEVALUATIONPOINT3D multiplies the coordinatevalues x, y, z in a
% suitable way to the matrix P (usualy the identity of suitable size)
% to create a highly efficient B-Spline evaluation
%
% Input:    Parametervalues of coordinates: x,y,z
%           Matrix:                         P (usually the identity)
%           Dimention of evaluation:        dim (3 for quadratic, 4 for
%                                           cubic)
% Output:   Matrix M, which can be multiplied to a vector of controlpoints
%           of a regular cell to get a volumetric evaluation
%
% computeEvaluationpoints3D(P,x,y,z,dim) ultiplies the coordinatevalues 
% x, y, z in a suitable way to the matrix P  to create a highly efficient 
% B-Spline evaluation


%compute necessary parameters
[~,b]=size(P);
[c,~]=size(x);
[e,~]=size(y);
[f,~]=size(z);

%compute matrix M
M=( ...
    reshape( ...
        permute( ...
            reshape( ...
                (z)*reshape( ...
                    permute( ...
                        permute( ...
                            reshape( ...
                                (y)*reshape( ...
                                    permute( ...
                                        reshape( ...
                                            (x)*reshape( ...
                                                P, ...
                                                [dim,dim*dim*b] ...
                                            ), ...
                                            [c,dim,dim,b] ...
                                        ), ...
                                        [2,1,3,4] ...
                                    ), ...
                                    [dim,c*dim*b] ...
                                ), ...
                                [e,c,dim,b] ...
                            ), ...
                            [2,1,3,4] ...
                        ), ...
                        [3,2,1,4] ...
                    ), ...
                    [dim,c*e*b] ...
                ), ...
                [f,e,c,b] ...
            ), ...
            [3,2,1,4] ...
        ), ...
        [c*e*f,b] ...
    ) ...
  );

end

