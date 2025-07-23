function [S,Lattice] = computeTriCubicSubdivisionMatrix(Input,varargin)
%COMPUTETRUCUBICSUBDIVISIONMATRIX computes out of a polytope structure the
%corresponding subdivision matrix, which refines a set of cubes sharing a
%central vertex to a smaller one.
%
% Input:    adjacency matrix or face matrix or edge matrix of the
%           corresponding dual polytope
%
% varargin: Several optional commands (see the manual)
%
% S:        Subdivision matrix
%
% Lattice:  Structure of the cubes. Each column represents the adjacency
%           matrix of one cube.
%
% For more information see the corresponding manual

%-------------------------------------------------------------------------
%       Input Parsing
%-------------------------------------------------------------------------


%The amount of additional input parameters  
NumberOfAdditionalInput = nargin - 1;

%Standard values for several options of varargin

status=false;
visualization=false;
PrintImages=false;
Prefix='';
Suffix='';
View=[0,20];

vararginDS=varargin;

%Checks, if the given keyword could be found in varargin. If yes, values
%were set and the keywords and the values were deleted from varargin. If
%they cannot be found, they are simply ignored.
i=1;
while i <= NumberOfAdditionalInput
    %Specifies whether a status of the algorithm progress should be shown
    if ischar(varargin{i}) && strcmp(varargin{i},'Status')
        status=true;
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;

    %Specifies whether a visualization of the algorithm progress should be shown
    elseif ischar(varargin{i}) && strcmp(varargin{i},'Visualization')
        visualization=true;
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;
    %Specifies if images should be saved
    elseif ischar(varargin{i}) && strcmp(varargin{i},'PrintImages')
        PrintImages=true;
        Prefix=varargin{i+1};
        Suffix=varargin{i+2};
        varargin(i)=[];
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-3;
    %Specifies angles of figures
    elseif ischar(varargin{i}) && strcmp(varargin{i},'View')
        View=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;

    else
        i=i+1;
    end
end

%-------------------------------------------------------------------------
%       Compute the 3-Polytope
%-------------------------------------------------------------------------


[PointCoordinates3D,KiteAll,StatusString,statusCounter] = Construct3Polytope(Input,vararginDS{:});

%Adding the central point
PointCoordinates3D=[PointCoordinates3D;[0,0,0]];

%The amount of control points
[AmountOfPoints,~]=size(PointCoordinates3D);

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') 3-Polytope was constructed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end



%-------------------------------------------------------------------------
%       Compute the Adjacency Matrix and the Face Matrix
%-------------------------------------------------------------------------


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

%The Face Matrix
Faces=computeFaceMatrix(AdjacencyMatrix,'false');

%The amount of faces
[AmountOfFaces,~]=size(Faces);

%The amount of cubes/volumes
AmountOfVolumes=length(AdjacencyMatrix);


%The amount of edges
AmountOfEdges=AmountOfPoints-AmountOfFaces-AmountOfVolumes-1;

%plots the structure
if visualization
    figure

    %get the amount of primal points
    AmountOfVertices=AmountOfVolumes;

    %project the dual ponts to the faces
    PointCoordinates3DTemp=PointCoordinates3D;
    PointCoordinates3DTemp(KiteAll(:,3),:)=1./(PointCoordinates3DTemp(KiteAll(:,3),1).^2+PointCoordinates3DTemp(KiteAll(:,3),2).^2+PointCoordinates3DTemp(KiteAll(:,3),3).^2).*PointCoordinates3DTemp(KiteAll(:,3),:);
    
    %plot all kites
    for i = 1:length(KiteAll)
           plot3(PointCoordinates3DTemp(KiteAll(i,[1,2]),1),PointCoordinates3DTemp(KiteAll(i,[1,2]),2),PointCoordinates3DTemp(KiteAll(i,[1,2]),3),'LineWidth',2,'Color',[0.5,0,0]);
           hold on
           plot3(PointCoordinates3DTemp(KiteAll(i,[2,3]),1),PointCoordinates3DTemp(KiteAll(i,[2,3]),2),PointCoordinates3DTemp(KiteAll(i,[2,3]),3),'LineWidth',2,'Color',[0,0,0.5]);
           plot3(PointCoordinates3DTemp(KiteAll(i,[3,4]),1),PointCoordinates3DTemp(KiteAll(i,[3,4]),2),PointCoordinates3DTemp(KiteAll(i,[3,4]),3),'LineWidth',2,'Color',[0,0,0.5]);
           plot3(PointCoordinates3DTemp(KiteAll(i,[4,1]),1),PointCoordinates3DTemp(KiteAll(i,[4,1]),2),PointCoordinates3DTemp(KiteAll(i,[4,1]),3),'LineWidth',2,'Color',[0.5,0,0]);
    end
    axis equal
    axis off

    %plot primal edge and dual points
    plot3(PointCoordinates3DTemp(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3DTemp(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3DTemp(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize',25);
    plot3(PointCoordinates3DTemp(1:AmountOfVertices,1),PointCoordinates3DTemp(1:AmountOfVertices,2),PointCoordinates3DTemp(1:AmountOfVertices,3),'.','Color',[0.5,0,0],'Markersize',25);
    plot3(PointCoordinates3DTemp(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,1),PointCoordinates3DTemp(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,2),PointCoordinates3DTemp(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,3),'.','Color',[0,0,0.5],'Markersize',25);


    %plot lines to central points
    for i = 1:length(KiteAll)
           plot3([PointCoordinates3DTemp(KiteAll(i,3),1),0],[PointCoordinates3DTemp(KiteAll(i,3),2),0],[PointCoordinates3DTemp(KiteAll(i,3),3),0],'LineWidth',2,'Color',[0.7,0.7,0]);     
    end

    %plot central point
    plot3(0,0,0,'.','Color',[0.5,0.5,0],'Markersize',30)

    title('Cubic Structure')
    
    view(View);
    
    %if the image should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"CubicStructure"+Suffix);
    end

end

%-------------------------------------------------------------------------
%       Compute the CDV Matrices of the Volumes (MBar)
%-------------------------------------------------------------------------

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Compute MBar.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%Projects the dual points to the primal faces
for i= AmountOfPoints-AmountOfFaces:AmountOfPoints-1
    PointCoordinates3D(i,:)=PointCoordinates3D(i,:)./(norm(PointCoordinates3D(i,:))^2);
end

%The number of the central point
CentralVertexNumber=length(PointCoordinates3D);

%The amount of edges, that lay on a primal vertex
VolumeGrade=zeros(AmountOfVolumes,1);
for i=1:AmountOfVolumes
    VolumeGrade(i)=sum(KiteAll(:,1)==i);
end

%The maximum amount of points in one volume
MaxPointsInVolume=2+2*max(VolumeGrade);

%The adjacency matrix for each volume
VolumeAdjacencyMatrix=zeros(MaxPointsInVolume,MaxPointsInVolume,AmountOfVolumes);

%The CDV matrix for each volume
VolumeColinDeVerdiereMatrix=zeros(MaxPointsInVolume,MaxPointsInVolume,AmountOfVolumes);

%Contains the information which volume contains which points
VolumePointIdentification=zeros(AmountOfVolumes,MaxPointsInVolume);

%Computes VolumeAdjacency, VolumeColinDeVerdiereMatrix and VolumeColinDeVerdiereMatrix
for i=1:AmountOfVolumes

    %Kites on the volume
    LocalKites=KiteAll(KiteAll(:,1)==i,:);

    %Set the primal point
    VolumePointIdentification(i,1)=i;

    %Set the edge points
    VolumePointIdentification(i,2*VolumeGrade(i)+2)=CentralVertexNumber;

    %Index of first edge point
    firstEdgePoint=2;
    EdgePointCounter=0;

    %Index of first face point
    firstFacePoint=VolumeGrade(i)+2;

    %index of central point
    CentralPoint=2*VolumeGrade(i)+2;

    %The amount of local kites
    [AmountOfLocalKites,~]=size(LocalKites);
    for j=1:AmountOfLocalKites

        %Goes through all kites and get an edge point
        PostitionEdgePoint1=find(VolumePointIdentification(i,:)==LocalKites(j,2));
        if isempty(PostitionEdgePoint1)
            VolumePointIdentification(i,firstEdgePoint+EdgePointCounter)=LocalKites(j,2);
            PostitionEdgePoint1=firstEdgePoint+EdgePointCounter;
            EdgePointCounter=EdgePointCounter+1;
        end

        %Goes through all kites and get a second edge point
        PostitionEdgePoint2=find(VolumePointIdentification(i,:)==LocalKites(j,4));
        if isempty(PostitionEdgePoint2)
            VolumePointIdentification(i,firstEdgePoint+EdgePointCounter)=LocalKites(j,4);
            PostitionEdgePoint2=firstEdgePoint+EdgePointCounter;
            EdgePointCounter=EdgePointCounter+1;
        end

        %Get the face point of the local kite
        VolumePointIdentification(i,firstFacePoint+j-1)=LocalKites(j,3);
        PostitionFacePoint=firstFacePoint+j-1;

        %Set the values of the adjacency matrix
        VolumeAdjacencyMatrix(1,PostitionEdgePoint1,i)=1;
        VolumeAdjacencyMatrix(PostitionEdgePoint1,1,i)=1;

        VolumeAdjacencyMatrix(1,PostitionEdgePoint2,i)=1;
        VolumeAdjacencyMatrix(PostitionEdgePoint2,1,i)=1;

        VolumeAdjacencyMatrix(PostitionFacePoint,PostitionEdgePoint1,i)=1;
        VolumeAdjacencyMatrix(PostitionEdgePoint1,PostitionFacePoint,i)=1;

        VolumeAdjacencyMatrix(PostitionFacePoint,PostitionEdgePoint2,i)=1;
        VolumeAdjacencyMatrix(PostitionEdgePoint2,PostitionFacePoint,i)=1;

        VolumeAdjacencyMatrix(PostitionFacePoint,CentralPoint,i)=1;
        VolumeAdjacencyMatrix(CentralPoint,PostitionFacePoint,i)=1;

    end

    %compute the volume CDV matrix
    LocalVertices=PointCoordinates3D(VolumePointIdentification(i,1:2*VolumeGrade(i)+2),:);
    LocalVertices=LocalVertices-1/2*LocalVertices(1,:);

    VolumeColinDeVerdiereMatrix(1:2*VolumeGrade(i)+2,1:2*VolumeGrade(i)+2,i)=computeColinDeVerdiereMatrix3D(LocalVertices,VolumeAdjacencyMatrix(1:2*VolumeGrade(i)+2,1:2*VolumeGrade(i)+2,i));
end

%Computes the matrix MBar
MBar=zeros(sum(2*VolumeGrade+2),AmountOfPoints);
LineCounter=1;

%Put the local matrices under each other to get MBar
for i=1:AmountOfVolumes
    MBar(LineCounter:LineCounter+2*VolumeGrade(i)+2-1,VolumePointIdentification(i,1:2*VolumeGrade(i)+2))=VolumeColinDeVerdiereMatrix(1:2*VolumeGrade(i)+2,1:2*VolumeGrade(i)+2,i);
    LineCounter=LineCounter+2*VolumeGrade(i)+2;
end

MBar=MBar./sum(MBar,2);
TildeMV=MBar*PointCoordinates3D;


%-------------------------------------------------------------------------
%       Compute the CDV Matrices of the Structure Elements (BBar)
%-------------------------------------------------------------------------

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Compute BBar.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%Initialise the matrix
BBar=zeros(AmountOfPoints,sum(2*VolumeGrade+2));

%More efficient to get the indices of the points
VolumePointIdentificationMBar=reshape(VolumePointIdentification',[AmountOfVolumes*MaxPointsInVolume,1]);
VolumePointIdentificationMBar(VolumePointIdentificationMBar==0)=[];

%set the entries for the points (just one entry is 1)
for i=1:AmountOfVolumes
    BBar(i,sum(2*VolumeGrade(1:i-1)+2)+1)=1;
end

%set the entries for the edges (just two entries are unequal to 0)
for i=AmountOfVolumes+1:AmountOfVolumes+AmountOfEdges

    %compute the two indices on the edge
    LinesOfTildeM=find(VolumePointIdentificationMBar==i);
    VolumeMidVertices=TildeMV(LinesOfTildeM,:);

    %get the values vor the CDV matrix
    alpha=abs(1/2*PointCoordinates3D(i,:)-VolumeMidVertices(2,:))/abs(VolumeMidVertices(1,:)-VolumeMidVertices(2,:));

    %set the values
    BBar(i,LinesOfTildeM(1))=alpha;
    BBar(i,LinesOfTildeM(2))=1-alpha;

end

%set the entries for the faces 
for i=AmountOfVolumes+AmountOfEdges+1:AmountOfVolumes+AmountOfEdges+AmountOfFaces

    %get the local face
    LocalFace=Faces(i-AmountOfVolumes-AmountOfEdges,:);
    LocalFace=LocalFace(LocalFace>0);
    
    %creates the adjacency matrix for the local face
    PolytopeAdjacency=create2PolytopeAdjacencyMatrix(length(LocalFace));
    Points=PointCoordinates3D(LocalFace,:);

    %gets the center of the face
    Center=PointCoordinates3D(i,:);
    
    %comutes the local CDV matrix
    [~,Weights]=computeCDVMatrix2D(PolytopeAdjacency,Points,Center);

    %adjust the weights
    Weights=sum(Weights);
    Weights=(Weights)/sum(Weights);

    %get the indices of the matrix
    LinesOfTildeM=VolumePointIdentificationMBar==i;
    
    %sort the weights
    [~,I]=sort(LocalFace);
    Weights=Weights(I);

    %set the entries of the matrix
    BBar(i,LinesOfTildeM)=Weights;
end


%set the entries for the volume

%get the indices
LinesOfTildeM=VolumePointIdentificationMBar==CentralVertexNumber;

%get the CDV matrix
CDV=computeColinDeVerdiereMatrix3D(PointCoordinates3D(1:AmountOfVolumes,:),AdjacencyMatrix);

%compute the weights
Weights=sum(CDV);
Weights=(Weights)/sum(Weights);

%set the entries
BBar(end,LinesOfTildeM)=Weights;

%-------------------------------------------------------------------------
%       Adjust the spectrum
%-------------------------------------------------------------------------

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Adjust the spectrum.\n'];
    fprintf(StatusString);
end

TildeC=BBar*MBar;
TildeC2=2*TildeC;


TildeR=zeros(length(TildeC));
TildeL=zeros(length(TildeC));


%Splits 2*TildeC into the sum of two matrices
for i=1:length(TildeC2)
    for j=1:length(TildeC2)
        if j>i
            TildeR(i,j)=TildeC2(i,j);
        elseif i>j
            TildeL(i,j)=TildeC2(i,j);
        end
    end
end

%split the diagonal entry
for i=1:length(TildeC2)
    TildeR(i,i)=1-sum(TildeR(i,:));
    TildeL(i,i)=1-sum(TildeL(i,:));
end

%compute the matrix exponential
MatrixR=expm((TildeR-eye(length(TildeR)))*log(2));
MatrixS=expm((TildeL-eye(length(TildeL)))*log(2));

%compute the subdivision matrix
S=MatrixS*MatrixR;

Lattice=computeCubicAdjacencyMatrix(KiteAll);

end



function [C] = computeColinDeVerdiereMatrix3D(Controlpoints,AdjacencyMatrix)
% COMPUTECOLINDEVERDIEREMATRIX3D comuputes the CDV matrix of the control
% points with given adjacency matrix
%
% Input:    Controlpoints: The control points that forms a 3-Polytope
%           AdjacencyMatrix: The geometric information, which point are
%           neighbored.
% Output:   C: CDV matrix
%
% computeColinDeVerdiereMatrix3D(Controlpoints,AdjacencyMatrix) computes the 
% CDV matrix of the control points with given adjacency matrix.

%Size of the adjacency matrix
[a,b]=size(AdjacencyMatrix);

%size of the control points
[c,d]=size(Controlpoints);

%The adjacency matrix has to be square
if a~=b
    error('A is not an adjacency matrix.');
end

%Adjacency matrix and controlpoints has to have the same size
if c~=a
    error('To few or less points in V compared with A.');
end

%Points has to be in R^3
if d~=3
    error('Points in V are not in R^3.');
end

%get the Face Matrix
F=computeFaceMatrix(AdjacencyMatrix,'false');

%compute the amount of faces
[AmountOfFaces,~]=size(F);

%The list of the dual points
DualPoints=zeros(AmountOfFaces,3);

%compute the dual points. The dual points has to lay in all dual faces,
%which are spanned by the vectors of the primal points
for i=1:AmountOfFaces
    n=nnz(F(i,:));
    GLS=Controlpoints(F(i,1:n),:);
    DualPoints(i,:)=(GLS'*GLS)\(GLS'*ones(n,1));
end

%the amount of primal points
[AmountOfVertices,~]=size(Controlpoints);

%The CDV matrix
C=zeros(AmountOfVertices);

%Graph of the adjacency matrix
G=graph(AdjacencyMatrix);

%The edges of the graph
PrimalEdges=G.Edges.EndNodes;

%the amount of edges
AmountOfEdges=length(PrimalEdges);

%Compute for every edge the entry of the CDV matrix
for i=1:AmountOfEdges

    %Get the point numbers for the primal points
    currentEdge=PrimalEdges(i,:);

    %Get the face, the two vertices have in common
    [FacesWithFirstVertex,~]=find(F==currentEdge(1));
    [FacesWithSecondVertex,~]=find(F==currentEdge(2));
    CommonFaces=FacesWithFirstVertex(ismember(FacesWithFirstVertex,FacesWithSecondVertex));
    
    %If there are more ore less then 2 faces somthing went wrong.
    if length(CommonFaces)~=2
        error('Something was wrong in the graph structure. Was the craph not planar or 3-connected?')
    end

    %The four points needed for calculation
    ui=Controlpoints(currentEdge(1),:);
    uj=Controlpoints(currentEdge(2),:);
    wf=DualPoints(CommonFaces(1),:);
    wg=DualPoints(CommonFaces(2),:);

    %Calculation of the matrix entry
    MatrixEntry=-abs(wf-wg)/abs(cross(ui,uj));

    %Setting the matrix entries
    C(currentEdge(1),currentEdge(2))=MatrixEntry;
    C(currentEdge(2),currentEdge(1))=MatrixEntry;
end

%Compute the sum of the off diagonal entries multiplied with the
%coordinates
SumMPointsM=C*Controlpoints;

%compute the coordinates of the diagonal entries
for i=1:AmountOfVertices
    C(i,i)=-SumMPointsM(i,:)/Controlpoints(i,:);
end


end


function [Lattice] = computeCubicAdjacencyMatrix(KiteAll)
%COMPUTECUBICADJACENCYMATRIX creates out of a 3-polytope a set of cubes.
%There is one cube for each corner point of the polytope.
%
%Input: KiteAll:    Information of the polytope in form of kites. This
%                   variable can created with the construct3Polytope
%                   function
%Output: Lattice:   The lattice consisting out of of the cubes. Each column
%                   is one adjacency matrix. It can be restored by a
%                   reshape command (the size is sqrt(Amount of Rows) x 
%                   sqrt(Amount of Rows)

%The amount of points and also the size of the adjacency matrix to compute
AmountOfPoints=max(max(KiteAll))+1;

%The amount of cubes
AmountOfLocalVolumes=length(unique(KiteAll(:,1)));

%creates one cube after another
for i=1:AmountOfLocalVolumes
   AdjacencyLocal=zeros(AmountOfPoints);

   %goes through each kite and looks if the corner point is part of the
   %kite
   for j = 1:length(KiteAll)

       %if so, the adjacency matrix ix computed
        if KiteAll(j,1)==i

            %computes the entries of the adjacency matrix
            AdjacencyLocal(i,KiteAll(j,2))=1;
            AdjacencyLocal(KiteAll(j,2),i)=1;

            AdjacencyLocal(i,KiteAll(j,4))=1;
            AdjacencyLocal(KiteAll(j,4),i)=1;


            AdjacencyLocal(KiteAll(j,3),KiteAll(j,2))=1;
            AdjacencyLocal(KiteAll(j,2),KiteAll(j,3))=1;

            AdjacencyLocal(KiteAll(j,3),KiteAll(j,4))=1;
            AdjacencyLocal(KiteAll(j,4),KiteAll(j,3))=1;


            AdjacencyLocal(KiteAll(j,3),end)=1;
            AdjacencyLocal(end,KiteAll(j,3))=1;
        end
   end

   %adds the adjacency matrix to the lattice
   if i==1
        Lattice=AdjacencyLocal;
   else
        Lattice=cat(3,Lattice,AdjacencyLocal);
   end
end


end






%Also a local function in computeBiCubicSubdivisionMatrix
function [CDVNormalized,CDV] = computeCDVMatrix2D(AdjacencyMatrix,Controlpoints,Center)
%COMPUTECDVMATRIX2D  computes a 2-dimensional CDV matrix for the given
%controlpints, whichs 2-polytope lies on a circle with mitpoint Center.
%
%Input:     AdjacencyMatrix:    The adjacency matrix of the structure
%           Controlpoints:      The controlpoints of the structure
%           Center:             The midpoint of the circle on which the
%                               2-polytope lies
% Output:   CDVNormalized:      The normalized CDV matrix
%           CDV:                The CDV matrix


%The amount of the adjacency matrix
[a,~]=size(AdjacencyMatrix);

%Initializing the CDV matrix
CDV=zeros(a);

%the size of the control points
[c,d]=size(Controlpoints);

%as this algorithm wors in R^3...
if d==2
    Controlpoints=[Controlpoints,zeros(c,1)];
end

%as this algorithm wors in R^3...
[c,d]=size(Center);
if d==2
    Center=[Center,zeros(c,1)];
end

%Shifts the control points such that the center is in [0,0,0]
PShifted=Controlpoints-Center;

%computes each entry
for i=1:a
    for j=1:a

        %entries are just unequal to 0, if the adjacency matrix is unequal
        %to 0
        if AdjacencyMatrix(i,j)==1

           %compute the entries
           Ui=[PShifted(i,:)];
           Uj=[PShifted(j,:)];
           UU=cross(Ui,Uj);
           CDV(i,j)=-1/(sqrt(UU*UU'));
  
        end
    end
    %unnecessary condition to avoid dividing by 0
    if sum(AdjacencyMatrix(i,:))>0

        %compute the duaginal entries
        UStrich=CDV(i,:)*PShifted;
        CDV(i,i)=-UStrich/PShifted(i,:);
    end

end

%normalize the CDV matrix
CDVNormalized=CDV./sum(CDV,2);


end







%Also a local function in computeBiCubicSubdivisionMatrix
function [AdjacencyMatrix] = create2PolytopeAdjacencyMatrix(n)
% CREATE2POLYTOPEADJACENCYMATRIX creates an adjacency matrix for a
% 2-polytope with n vertices
%
% Input:    n:                  amount of vertices
% Output:   AdjacencyMatrix:    the adjacency matrix
%
% create2PolytopeAdjacencyMatrix(n) creates an adjacency matrix for a
% 2-polytope with n vertices

%initializing the adjacency matrix
AdjacencyMatrix=zeros(n);

%goes through each entry
for i=1:n
    for j=1:n

        %if one of the two conditions is fulfilled the entry is 1
        if j==i+1
            AdjacencyMatrix(i,j)=1;
        end
        if j==i-1
            AdjacencyMatrix(i,j)=1;
        end
    end
end

%the two special cases (mod n)
AdjacencyMatrix(1,n)=1;
AdjacencyMatrix(n,1)=1;
end