function [S,Lattice] = computeTriQuadraticSubdivisionMatrixV2Big(Input,varargin)
%COMPUTETRUQUADRATICSUBDIVISIONMATRIXV2BIG computes out of a polytope structure the
%corresponding subdivision matrix, which refines a shell to a smaller one.
%
% Input:    adjacency matrix or face matrix or edge matrix of the
%           polytope
%
% varargin: Several optional commands (see the manual)
%
% S:        Subdivision matrix
%
% Lattice:  Structure of the polytopes. Each column represents the adjacency
%           matrix of one polytoe.
%
% For more information see the corresponding manual


%-------------------------------------------------------------------------
%       Input Parsing
%-------------------------------------------------------------------------

%The amount of additional input parameters  
NumberOfAdditionalInput = nargin - 1;

%The additional options for the central computation of the subdivision
%algorithm
vararginSmallCentral={};

%The additional options for the computation of all other parts of the
%subdivision algorithm. Input check is here not neccecary as all outer
%elements are prisms by construction.
vararginSmallOutside={'PreventInputCheck'};

status=false;
StatusString="";
statusCounter=1;
KeepFaceOrder=false;
visualization=false;
PrintImages=false;
Prefix='';
Suffix='';
View=[0,20];

%Checks, if the given keyword could be found in varargin. If yes, values
%were set and the keywords and the values were deleted from varargin. If
%they cannot be found, they are simply ignored.
i=1;
while i <= NumberOfAdditionalInput

    %Specifies whether the eigenstructure shult be plottet or not
    if ischar(varargin{i}) && strcmp(varargin{i},'Visualization')
        visualization=true;
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;

    %Specifies if images should be saved
    elseif ischar(varargin{i}) && strcmp(varargin{i},'PrintImages')
        PrintImages=true;
        Prefix=varargin{i+1};
        Suffix=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;
    %Specifies angles of figures
    elseif ischar(varargin{i}) && strcmp(varargin{i},'View')
        View=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;

    

    %Specifies whether a status of the algorithm progress should be shown
    elseif ischar(varargin{i}) && strcmp(varargin{i},'Status')
        status=true;
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;

    %Specifies if the central subdivision algorithm should do an input
    %check
    elseif ischar(varargin{i}) && strcmp(varargin{i},'PreventInputCheck')
        vararginSmallCentral=[vararginSmallCentral(:)',{'PreventInputCheck'}];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;

    %Specifies the tolerance of the numerical parts of the algorithm
    elseif ischar(varargin{i}) && strcmp(varargin{i},'Tolerance')
        vararginSmallCentral=[vararginSmallCentral(:)',{'Tolerance'},varargin{i+1}];
        vararginSmallOutside=[vararginSmallOutside(:)',{'Tolerance'},varargin{i+1}];
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;

    %Specifies the amount of iterations of the numerical parts of the
    %algorithm
    elseif ischar(varargin{i}) && strcmp(varargin{i},'MaxIterations')
        if varargin{i+1} < 1
            vararginSmallCentral=[vararginSmallCentral(:)',{'MaxIterations'},100];
            vararginSmallOutside=[vararginSmallOutside(:)',{'MaxIterations'},100];
            warning('MaxIterations has to be a positive number. Set the value to 100.')
        else
            vararginSmallCentral=[vararginSmallCentral(:)',{'MaxIterations'},varargin{i+1}];
            vararginSmallOutside=[vararginSmallOutside(:)',{'MaxIterations'},varargin{i+1}];
        end
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;

    %Specifies if the order of the input faces should kept (only if Input
    %is a face matrix)
    elseif ischar(varargin{i}) && strcmp(varargin{i},'KeepFaceOrder')
        KeepFaceOrder=true;
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;
    else
        i=i+1;
    end
end

%-------------------------------------------------------------------------
%       Initialization
%-------------------------------------------------------------------------

%Checks whether the input is a face matrix, an adjacency matrix or an edge matrix 
%also check some basic plausibilitys if the input is valid
%creates then the adjacency matrix for all Input types

[rows,cols]=size(Input);

%If the input has just 2 columns then we have an edge matrix
if cols ==2
    InputType='Edge';
    EdgeMatrix=Input;

    %The amount of vertices
    AmountOfVertices=max(max(EdgeMatrix));
    
    %Comparing the labeling of the vertices with 1,2,...,AmountOfVertices
    if sum((1:AmountOfVertices)==(unique(EdgeMatrix))') ~= AmountOfVertices
        error('Vertices are not named uniformly.');
    end

    G=graph(EdgeMatrix(:,1),EdgeMatrix(:,2));
    AdjacencyMatrixSmall=full(G.adjacency);

    
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Edge-List input type was identified.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end

%If the input is square and every entry is either 0 or 1, then we have an adjacency matrix 
elseif rows==cols && sum(sum((Input==0)+(Input==1)))==rows^2
    InputType='Adjacency';
    AdjacencyMatrixSmall=Input;
    

    if sum(sum(AdjacencyMatrixSmall==AdjacencyMatrixSmall'))~=rows^2
        error('Adjacency matrix as input must be symmetric.')
    end
    
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Adjacency matrix input type was identified.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end

%Otherwise we have a face matrix (or something invalid)
else
    InputType='Face';

    %The amount of vertices
    AmountOfVertices=max(max(Input));
    
    %Comparing the labeling of the vertices with 1,2,...,AmountOfVertices
    if sum((1:AmountOfVertices)==(unique(Input(Input>0)))')~=AmountOfVertices
        error('Vertices are not named uniformly.');
    end

    
    AdjacencyMatrixSmall=FaceToAdjacencyMatrix(Input);

    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Face matrix input type was identified.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end
end




%The orientation of the blocks depend on the numeration of the faces. If
%the input was a face matrix, the user can determine, if it want to keep
%this order, of the the algorithm should sort it on its own

%Take algorithm order
if ~KeepFaceOrder
    
    %List of faces of the polytope
    FaceMatrix=computeFaceMatrix(AdjacencyMatrixSmall,false);

%Want manual order but not set it. Take algorithm order
elseif KeepFaceOrder && ~strcmp(InputType,'Face')
    warning('Face order should be kept, but input was not a face matrix. Take new face order.')

    %List of Faces of the Polytope
    FaceMatrix=computeFaceMatrix(AdjacencyMatrixSmall);

%Take manual order
else
    %List of Faces of the Polytope
    FaceMatrix=Input;
end

G=graph(AdjacencyMatrixSmall);

%A List of all edges of the adjacency matrix.
EdgeList=G.Edges.EndNodes;

%The amount of blocks is the amount of vertices in the inner polygon
%(for each vertex there is one block).
AmountOfBlocks=max(max(FaceMatrix));

%Amount of faces of the polygon
[AmountOfFaces,LargestFace]=size(FaceMatrix);

%The amount of edges of the polygon 
AmountOfEdges=sum(sum(AdjacencyMatrixSmall))/2;

%The size of the block depend on the shape of the trapezohedron. A
%trapezohedron has 2n+2 vertices, 4n edges, 2n faces and 1 volume. Due to
%refinement a block has 8n+3 vertices (for a cube that are 27)
BlockSize = zeros(AmountOfBlocks,1);

%The amount of faces on which the point of the block lay in the small
%volume.
BlockFaces= zeros(AmountOfBlocks,1);

%The amount of volumes for each block in a lattice. This is equal to the
%amount of points in the corresponding trapezohedron, so it is 2n+2. For a
%cube, it is 8.
BlockVolumes = zeros(AmountOfBlocks,1);

%Compute the variables for each block
for i = 1:AmountOfBlocks
    faces=find(FaceMatrix==i);
    BlockFaces(i)=length(faces);
    BlockSize(i)=length(faces)*8+3;
    BlockVolumes(i)=length(faces)*2+2;
end

%the FaceBlockList has the information, which block lays on which faces. As
%the sequence is important, it is not that easy to compute
FaceBlockList=zeros(AmountOfBlocks,max(BlockFaces));

%compute the FaceBlockLost
for i=1:AmountOfBlocks

    %starts for each Block with the minimum face
    [row,cols]=find(FaceMatrix==i);
    
    [nextFace,position]=min(row);
    ColOfNextFace=cols(position);
    LastAdditionalEntry=0;
    
    FaceCounter=1;
    while nextFace ~= FaceBlockList(i,1)
        FaceBlockList(i,FaceCounter)=nextFace;
        actualFace=nextFace;
        actualCol=ColOfNextFace;
        
        %Determines, if there is a col break (last entry and first entry of
        %a line are neighbored)
        if actualCol == LargestFace || FaceMatrix(actualFace,actualCol+1) == 0

            %checks which line in the face was seen in the last loop. Take
            %the other one
            if FaceMatrix(actualFace,actualCol-1)== LastAdditionalEntry || LastAdditionalEntry == 0
                AdditionalEntry= FaceMatrix(actualFace,1);
            else
                AdditionalEntry= FaceMatrix(actualFace,actualCol-1);
            end

        %Determines, if there is a col break (last entry and first entry of
        %a line are neighbored)
        elseif actualCol==1

            %checks which line in the face was seen in the last loop. Take
            %the other one
            AmountOfEntries=nnz(FaceMatrix(actualFace,:));
            if FaceMatrix(actualFace,AmountOfEntries)== LastAdditionalEntry || LastAdditionalEntry == 0
                AdditionalEntry= FaceMatrix(actualFace,actualCol+1);
            else
                AdditionalEntry= FaceMatrix(actualFace,AmountOfEntries);
            end

        %if there is no col break
        else

            %checks which line in the face was seen in the last loop. Take
            %the other one
            if FaceMatrix(actualFace,actualCol-1)== LastAdditionalEntry || LastAdditionalEntry == 0
                AdditionalEntry= FaceMatrix(actualFace,actualCol+1);
            else
                AdditionalEntry= FaceMatrix(actualFace,actualCol-1);
            end
        end
        
        %Searches for the faces, which have also the aditional entry
        [rowAdditional,~]=find(FaceMatrix==AdditionalEntry);

        %as the current face also have both entry one has to dertermine,
        %which one are the two with both entries
        SameFaces=ismember(row,rowAdditional);
        SameFaces=row(SameFaces);
    
        %searches for the next face (the face, which has both entries and
        %is not the current face=
        if SameFaces(1)== actualFace
            nextFace=SameFaces(2);
        else
            nextFace=SameFaces(1);
        end

        position=row==nextFace;
        ColOfNextFace=cols(position);
        LastAdditionalEntry=AdditionalEntry;
    
        FaceCounter=FaceCounter+1;

    end
   
end

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') FaceBlockList was computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end


%The amount of lattice volumes (one for the initial volume, 2 for each face, 4 for each edge
% (in the initial volume) and the amount of volumes per block)
AmountOfVolumesBig=1+2*AmountOfFaces+4*AmountOfEdges+sum(BlockVolumes);


%Initial values for the Lattice and the Subdivision matrix. As both are
%sparse, we save just row, col and value
LatticeCounter=1;
Indexcounter=0;

LatticeRow=zeros(length(AdjacencyMatrixSmall)^2+100000,1);
LatticeCol=zeros(length(AdjacencyMatrixSmall)^2+100000,1);
LatticeValue=zeros(length(AdjacencyMatrixSmall)^2+100000,1);

SBigRow=zeros(length(AdjacencyMatrixSmall)^2+100000,1);
SBigCol=zeros(length(AdjacencyMatrixSmall)^2+100000,1);
SBigValue=zeros(length(AdjacencyMatrixSmall)^2+100000,1);

%-------------------------------------------------------------------------
%       Center Subdivision
%-------------------------------------------------------------------------

%Compute the subdivision weights of the central element
SCenter=computeTriQuadraticSubdivisionMatrixV2(FaceMatrix,vararginSmallCentral{:});


%compute the Index of the entries
I2=(1:AmountOfBlocks)';
I4=(1:AmountOfBlocks)';

I2=repmat(I2,AmountOfBlocks,1);
I4=repelem(I4,AmountOfBlocks,1);
I1=ones(length(I2),1);
I3=ones(length(I2),1);

P=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));


%prelocating more space
if length(P)+Indexcounter>length(LatticeRow)
    [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
end

%Set the computed values
LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
LatticeValue(Indexcounter+1:Indexcounter+length(P))=reshape(AdjacencyMatrixSmall,[length(P),1]);

SBigRow(Indexcounter+1:Indexcounter+length(P))=P;
SBigCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1);
SBigValue(Indexcounter+1:Indexcounter+length(P))=reshape(SCenter,[length(P),1]);

Indexcounter=Indexcounter+length(P);

LatticeCounter=LatticeCounter+1;

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Center subdivision was computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end


%-------------------------------------------------------------------------
%       Face Subdivision
%-------------------------------------------------------------------------

%compute the maximum face size
maxN=0;
for i=1:AmountOfFaces
    n=nnz(FaceMatrix(i,:));
    if n>maxN
        maxN=n;
    end
end

PreconditionM=zeros(2*maxN,2*maxN,maxN);
seen=zeros(1,maxN);

%computes the face matrix for each prims valence, that is needed to improve
%speed
for i=1:AmountOfFaces
    n=nnz(FaceMatrix(i,:));
    if seen(n)==0
        PrismaPolygon=computePrismFaceMatrix(n);
        M=computeTriQuadraticSubdivisionMatrixV2(PrismaPolygon,'PreventInputCheck');
        PreconditionM(1:2*n,1:2*n,n)=M;
        seen(n)=1;
    end
end


%goes through all faces and compute the values of the matrix
for i=1:AmountOfFaces

    %initializing local varibales
    P1Temp=zeros(100000,1);
    P2Temp=zeros(100000,1);
    PValue1Temp=zeros(100000,1);
    PValue2Temp=zeros(100000,1);
    SBigTemp=zeros(100000,1);
    SBigValueTemp=zeros(100000,1);
    localCounter=0;
    
    %valence of current face
    n=nnz(FaceMatrix(i,:));

    %faceMatrix of current face
    PrismaPolygon=computePrismFaceMatrix(n);

    %current local subdivision matrix
    M=PreconditionM(1:2*n,1:2*n,n);

    %current local adjacency matrix matrix
    MAdjacency=FaceToAdjacencyMatrix(PrismaPolygon);

    %goes for each block on the face through each row
    for j=1:n

        %set the current row
        currentBlockRow=FaceMatrix(i,j);
        FaceNumberRow=find(FaceBlockList(currentBlockRow,:)==i);
        
        %set indicex
        EntriesRow=[1+2*(FaceNumberRow-1)+1,1+2*(BlockFaces(currentBlockRow))+1+4*(FaceNumberRow-1)+1];
        EntriesRowJ=[1+2*(FaceNumberRow-1)+1,1+2*(BlockFaces(currentBlockRow))+1+4*(FaceNumberRow-1)+1];
        EntriesColJ=[1, 1+2*(FaceNumberRow-1)+1];

        %goes for each block on the face through each col
        for k=1:n

            %set the current col
            currentBlockCol=FaceMatrix(i,k);
            FaceNumberCol=find(FaceBlockList(currentBlockCol,:)==i);

            %set the indices
            EntriesRowK=[1+2*(FaceNumberCol-1)+1,1+2*(BlockFaces(currentBlockCol))+1+4*(FaceNumberCol-1)+1];
            EntriesCol=[1, 1+2*(FaceNumberCol-1)+1];
            EntriesColK=[1, 1+2*(FaceNumberCol-1)+1];
            I1=EntriesRow';
            I3=EntriesCol';
            I1=repmat(I1,2,1);
            I3=repelem(I3,2,1);
            I2=ones(length(I1),1)*currentBlockRow;
            I4=ones(length(I1),1)*currentBlockCol;
            NewIndizesS= computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));

            I1=EntriesRowJ';
            I3=EntriesRowK';
            I1=repmat(I1,2,1);
            I3=repelem(I3,2,1);
            I2=ones(length(I1),1)*currentBlockRow;
            I4=ones(length(I1),1)*currentBlockCol;
            NewIndizesP1=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));

            I1=EntriesColJ';
            I3=EntriesColK';
            I1=repmat(I1,2,1);
            I3=repelem(I3,2,1);
            I2=ones(length(I1),1)*currentBlockRow;
            I4=ones(length(I1),1)*currentBlockCol;
            NewIndizesP2=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));

            %allocate space
            maxNewEntries=max([length(NewIndizesS),length(NewIndizesP1),length(NewIndizesP2)]);
            if maxNewEntries+localCounter>length(SBigTemp)
                SBigTempT=SBigTemp;
                SBigValueTempT=SBigValueTemp;
                P1TempT=P1Temp;
                PValue1TempT=PValue1Temp;
                P2TempT=P2Temp;
                PValue2TempT=PValue2Temp;

                SBigTemp=zeros(length(SBigTemp)+100000+maxNewEntries,1);
                SBigValueTemp=zeros(length(SBigValueTemp)+100000+maxNewEntries,1);
                P1Temp=zeros(length(P1Temp)+100000+maxNewEntries,1);
                PValue1Temp=zeros(length(PValue1Temp)+100000+maxNewEntries,1);
                P2Temp=zeros(length(P2Temp)+100000+maxNewEntries,1);
                PValue2Temp=zeros(length(PValue2Temp)+100000+maxNewEntries,1);

                SBigTemp(1:length(SBigTempT))=SBigTempT;
                SBigValueTemp(1:length(SBigTempT))=SBigValueTempT;
                P1Temp(1:length(SBigTempT))=P1TempT;
                PValue1Temp(1:length(SBigTempT))=PValue1TempT;
                P2Temp(1:length(SBigTempT))=P2TempT;
                PValue2Temp(1:length(SBigTempT))=PValue2TempT;

            end

            %set the values
            SBigTemp(localCounter+1:localCounter+length(NewIndizesS))=NewIndizesS;
            SBigValueTemp(localCounter+1:localCounter+length(NewIndizesS))=reshape(M([j,j+n],[k,k+n]),[4,1]);

            P1Temp(localCounter+1:localCounter+length(NewIndizesP1))=NewIndizesP1;
            PValue1Temp(localCounter+1:localCounter+length(NewIndizesP1))=reshape(MAdjacency([j,j+n],[k,k+n]),[4,1]);

            P2Temp(localCounter+1:localCounter+length(NewIndizesP2))=NewIndizesP2;
            PValue2Temp(localCounter+1:localCounter+length(NewIndizesP2))=reshape(MAdjacency([j,j+n],[k,k+n]),[4,1]);
            
            localCounter=localCounter+maxNewEntries;
            
            

            
        end
    end
    
    %reduce variable size
    SBigTemp=SBigTemp(1:localCounter);
    SBigValueTemp=SBigValueTemp(1:localCounter);
    P1Temp=P1Temp(1:localCounter);
    PValue1Temp=PValue1Temp(1:localCounter);
    P2Temp=P2Temp(1:localCounter);
    PValue2Temp=PValue2Temp(1:localCounter);

    %set the values to the big matrices
    P=P1Temp;

    %allocate space
    if max(length(P),length(SBigTemp))+Indexcounter>length(LatticeRow)
        if length(P)>length(SBigTemp)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        else
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,SBigTemp);
        end
    end

    %set the values
    LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
    LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
    LatticeValue(Indexcounter+1:Indexcounter+length(P))=PValue1Temp;
    
    SBigRow(Indexcounter+1:Indexcounter+length(SBigTemp))=SBigTemp;
    SBigCol(Indexcounter+1:Indexcounter+length(SBigTemp))=ones(length(SBigTemp),1);
    SBigValue(Indexcounter+1:Indexcounter+length(SBigTemp))=SBigValueTemp;
    
    Indexcounter=Indexcounter+max(length(P),length(SBigTemp));
    LatticeCounter=LatticeCounter+1;

    %set the values to the big matrices
    P=P2Temp;

    %allocate space
    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end

    %set the values
    LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
    LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
    LatticeValue(Indexcounter+1:Indexcounter+length(P))=PValue2Temp;
    
    Indexcounter=Indexcounter+length(P);
    LatticeCounter=LatticeCounter+1;
end


if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Face subdivision was computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%-------------------------------------------------------------------------
%       Edge Subdivision
%-------------------------------------------------------------------------

%compute the local subdivision and adjacency matrix
PrismaPolygon=computePrismFaceMatrix(4);
M=computeTriQuadraticSubdivisionMatrixV2(PrismaPolygon,'PreventInputCheck');
MAdjacency=FaceToAdjacencyMatrix(PrismaPolygon);

%goes through all edges and compute the entries
for i=1:AmountOfEdges

    %the current edge
    currentEdge=EdgeList(i,:);

    %the two faces on this edge
    SameFaces=FaceBlockList(currentEdge(1),ismember(FaceBlockList(currentEdge(1),:), FaceBlockList(currentEdge(2),:)));
    Face1=SameFaces(1);
    Face2=SameFaces(2);
    n=nnz(FaceBlockList(currentEdge(1),:));
    Position11= find(FaceBlockList(currentEdge(1),:)==Face1);
    Position12= find(FaceBlockList(currentEdge(1),:)==Face2);
    Swap1=false;

    %depending on the index some values must be swapped
    if Position11 == 1 && Position12 == n
        Position11=Position12;
        Swap1=true;
    elseif Position11 == n && Position12 == 1

    elseif Position11>Position12
        Position11=Position12;
        Swap1=true;
    end
    
    % if there was a swap
    if Swap1

        %compute the indices
        EntriesRow1=[1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+4,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+3,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+2,...
                    1+2*(Position11-1)+2];

        EntriesLeft1=[1+2*(Position11-1)+2,...
                1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+2,...
                1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+1,...
                1+2*(Position11-1)+1];

        
        if Position11 == n
            EntriesCol1=[2,...
                        1+2*(Position11-1)+2,...
                        1+2*(Position11-1)+1,...
                        1];

            EntriesRight1=[2,...
                    1+2*(BlockFaces(currentEdge(1)))+1+1,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+4,...
                    1+2*(Position11-1)+2];
        else
    
            EntriesCol1=[1+2*(Position11-1)+3,...
                        1+2*(Position11-1)+2,...
                        1+2*(Position11-1)+1,...
                        1];

            EntriesRight1=[1+2*(Position11-1)+3,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+5,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+4,...
                    1+2*(Position11-1)+2];
        end


    else

        %compute the indices
        EntriesRow1=[1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+2,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+3,...
                    1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+4,...
                    1+2*(Position11-1)+2];

        EntriesLeft1=[1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+1,...
                1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+2,...
                1+2*(Position11-1)+2,...
                1+2*(Position11-1)+1];
    
        if Position11 == n
            EntriesCol1=[1+2*(Position11-1)+1,...
                        1+2*(Position11-1)+2,...
                        2,...
                        1];

            EntriesRight1=[1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+4,...
                1+2*(BlockFaces(currentEdge(1)))+1+1,...
                2,...
                1+2*(Position11-1)+2];
        else
    
            EntriesCol1=[1+2*(Position11-1)+1,...
                        1+2*(Position11-1)+2,...
                        1+2*(Position11-1)+3,...
                        1];

            EntriesRight1=[1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+4,...
                1+2*(BlockFaces(currentEdge(1)))+1+4*(Position11-1)+5,...
                1+2*(Position11-1)+3,...
                1+2*(Position11-1)+2];
        end

    end

    %do the same for the second one
    n=nnz(FaceBlockList(currentEdge(2),:));


    Position21= find(FaceBlockList(currentEdge(2),:)==Face1);
    Position22= find(FaceBlockList(currentEdge(2),:)==Face2);
    
    %for some indices a swap is needed
    Swap2=false;
    if Position21 == 1 && Position22 == n
        Position21=Position22;
        Swap2=true;
    elseif Position21 == n && Position22 == 1

    elseif Position21>Position22
        Position21=Position22;
        Swap2=true;
    end

    %if a swap was needed
    if Swap2

        %compute the indices
        EntriesRow2=[1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+4,...
                    1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+3,...
                    1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+2,...
                    1+2*(Position21-1)+2];

        EntriesLeft2=[1+2*(Position21-1)+2,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+2,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+1,...
                1+2*(Position21-1)+1];
    
        if Position21 == n
            EntriesCol2=[2,...
                    1+2*(Position21-1)+2,...
                    1+2*(Position21-1)+1,...
                    1];

            EntriesRight2=[2,...
                1+2*(BlockFaces(currentEdge(2)))+1+1,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+4,...
                1+2*(Position21-1)+2];
        else
    
            EntriesCol2=[1+2*(Position21-1)+3,...
                        1+2*(Position21-1)+2,...
                        1+2*(Position21-1)+1,...
                        1];

            EntriesRight2=[1+2*(Position21-1)+3,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+5,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+4,...
                1+2*(Position21-1)+2];
        end
    else
    
        %compute the indices
        EntriesRow2=[1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+2,...
                    1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+3,...
                    1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+4,...
                    1+2*(Position21-1)+2];

        EntriesLeft2=[1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+1,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+2,...
                1+2*(Position21-1)+2,...
                1+2*(Position21-1)+1];
    
        if Position21 == n
            EntriesCol2=[1+2*(Position21-1)+1,...
                    1+2*(Position21-1)+2,...
                    2,...
                    1];

            EntriesRight2=[1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+4,...
                1+2*(BlockFaces(currentEdge(2)))+1+1,...
                2,...
                1+2*(Position21-1)+2];

        else
    
            EntriesCol2=[1+2*(Position21-1)+1,...
                        1+2*(Position21-1)+2,...
                        1+2*(Position21-1)+3,...
                        1];

            EntriesRight2=[1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+4,...
                1+2*(BlockFaces(currentEdge(2)))+1+4*(Position21-1)+5,...
                1+2*(Position21-1)+3,...
                1+2*(Position21-1)+2];
        end

    end
           

    %Set the indieces for the big matrix
    I1=EntriesRow1';
    I3=EntriesCol1';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(1);
    I4=ones(length(I1),1)*currentEdge(1);
    P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value1=reshape(M(1:4,1:4),[16,1]);

    I1=EntriesRow1';
    I3=EntriesCol2';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(1);
    I4=ones(length(I1),1)*currentEdge(2);
    P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value2=reshape(M(1:4,5:8),[16,1]);

    I1=EntriesRow2';
    I3=EntriesCol1';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(2);
    I4=ones(length(I1),1)*currentEdge(1);
    P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value3=reshape(M(5:8,1:4),[16,1]);

    I1=EntriesRow2';
    I3=EntriesCol2';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(2);
    I4=ones(length(I1),1)*currentEdge(2);
    P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value4=reshape(M(5:8,5:8),[16,1]);


    P=[P21;P22;P23;P24];
    P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];

    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end
    
   
    %set the indices
    SBigRow(Indexcounter+1:Indexcounter+length(P))=P;
    SBigCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1);
    SBigValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
    
    Indexcounter=Indexcounter+length(P);


    %Set the indieces for the big matrix
    I1=EntriesRow1';
    I3=EntriesRow1';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(1);
    I4=ones(length(I1),1)*currentEdge(1);
    P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value1=reshape(MAdjacency(1:4,1:4),[16,1]);

    I1=EntriesRow1';
    I3=EntriesRow2';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(1);
    I4=ones(length(I1),1)*currentEdge(2);
    P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value2=reshape(MAdjacency(1:4,5:8),[16,1]);

    I1=EntriesRow2';
    I3=EntriesRow1';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(2);
    I4=ones(length(I1),1)*currentEdge(1);
    P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value3=reshape(MAdjacency(5:8,1:4),[16,1]);

    I1=EntriesRow2';
    I3=EntriesRow2';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(2);
    I4=ones(length(I1),1)*currentEdge(2);
    P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value4=reshape(MAdjacency(5:8,5:8),[16,1]);
    
    P=[P21;P22;P23;P24];
    P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];

    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end
    
    %set the indices
    LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
    LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
    LatticeValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
    

    Indexcounter=Indexcounter+length(P);
    LatticeCounter=LatticeCounter+1;

  
    %Set the indices for the big matrix
    I1=EntriesCol1';
    I3=EntriesCol1';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(1);
    I4=ones(length(I1),1)*currentEdge(1);
    P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value1=reshape(MAdjacency(1:4,1:4),[16,1]);

    I1=EntriesCol1';
    I3=EntriesCol2';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(1);
    I4=ones(length(I1),1)*currentEdge(2);
    P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value2=reshape(MAdjacency(1:4,5:8),[16,1]);

    I1=EntriesCol2';
    I3=EntriesCol1';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(2);
    I4=ones(length(I1),1)*currentEdge(1);
    P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value3=reshape(MAdjacency(5:8,1:4),[16,1]);

    I1=EntriesCol2';
    I3=EntriesCol2';
    I1=repmat(I1,4,1);
    I3=repelem(I3,4,1);
    I2=ones(length(I1),1)*currentEdge(2);
    I4=ones(length(I1),1)*currentEdge(2);
    P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P2Value4=reshape(MAdjacency(5:8,5:8),[16,1]);




    P=[P21;P22;P23;P24];
    P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];

    
    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end

    %set the indices
    LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
    LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
    LatticeValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
    

    Indexcounter=Indexcounter+length(P);
    LatticeCounter=LatticeCounter+1;

    %The next values depends, if there was two or no swap
    if (Swap1 && Swap2) || (~Swap1 && ~ Swap2)
    
        %Set the indices for the big matrix
        I1=EntriesLeft1';
        I3=EntriesLeft1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(1);
        P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value1=reshape(MAdjacency(1:4,1:4),[16,1]);
    
        I1=EntriesLeft1';
        I3=EntriesLeft2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(2);
        P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value2=reshape(MAdjacency(1:4,5:8),[16,1]);
    
        I1=EntriesLeft2';
        I3=EntriesLeft1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(1);
        P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value3=reshape(MAdjacency(5:8,1:4),[16,1]);
    
        I1=EntriesLeft2';
        I3=EntriesLeft2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(2);
        P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value4=reshape(MAdjacency(5:8,5:8),[16,1]);
    
    
        P=[P21;P22;P23;P24];
        P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];
    
    
        if length(P)+Indexcounter>length(LatticeRow)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        end
        
        %set the indices
        LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
        LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
        LatticeValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
        
        
        Indexcounter=Indexcounter+length(P);
        LatticeCounter=LatticeCounter+1;
    
    
        %Set the indices for the big matrix
        I1=EntriesRight1';
        I3=EntriesRight1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(1);
        P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value1=reshape(MAdjacency(1:4,1:4),[16,1]);
    
        I1=EntriesRight1';
        I3=EntriesRight2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(2);
        P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value2=reshape(MAdjacency(1:4,5:8),[16,1]);
    
        I1=EntriesRight2';
        I3=EntriesRight1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(1);
        P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value3=reshape(MAdjacency(5:8,1:4),[16,1]);
    
        I1=EntriesRight2';
        I3=EntriesRight2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(2);
        P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value4=reshape(MAdjacency(5:8,5:8),[16,1]);
    
    
    
    
        P=[P21;P22;P23;P24];
        P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];
    
        if length(P)+Indexcounter>length(LatticeRow)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        end
        
        %set the indices
        LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
        LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
        LatticeValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
        
       
        
        Indexcounter=Indexcounter+length(P);
        LatticeCounter=LatticeCounter+1;
    
    %if there was only one swap
    else

        if Swap2
            Temp=EntriesLeft1(1);
            EntriesLeft1(1)=EntriesLeft1(3);
            EntriesLeft1(3)=Temp;
    
            Temp=EntriesLeft1(1);
            EntriesLeft1(1)=EntriesLeft1(4);
            EntriesLeft1(4)=Temp;
    
            Temp=EntriesLeft1(2);
            EntriesLeft1(2)=EntriesLeft1(3);
            EntriesLeft1(3)=Temp;
    
            Temp=EntriesRight1(1);
            EntriesRight1(1)=EntriesRight1(3);
            EntriesRight1(3)=Temp;
    
            Temp=EntriesRight1(1);
            EntriesRight1(1)=EntriesRight1(4);
            EntriesRight1(4)=Temp;
    
            Temp=EntriesRight1(2);
            EntriesRight1(2)=EntriesRight1(3);
            EntriesRight1(3)=Temp;

        else
            Temp=EntriesLeft1(1);
            EntriesLeft1(1)=EntriesLeft1(3);
            EntriesLeft1(3)=Temp;

            Temp=EntriesLeft1(1);
            EntriesLeft1(1)=EntriesLeft1(2);
            EntriesLeft1(2)=Temp;
    
            Temp=EntriesLeft1(4);
            EntriesLeft1(4)=EntriesLeft1(3);
            EntriesLeft1(3)=Temp;
    
            Temp=EntriesRight1(1);
            EntriesRight1(1)=EntriesRight1(3);
            EntriesRight1(3)=Temp;

            Temp=EntriesRight1(1);
            EntriesRight1(1)=EntriesRight1(2);
            EntriesRight1(2)=Temp;
    
            Temp=EntriesRight1(4);
            EntriesRight1(4)=EntriesRight1(3);
            EntriesRight1(3)=Temp;
    
           
        end

        %Set the indices for the big matrix
        I1=EntriesLeft1';
        I3=EntriesLeft1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(1);
        P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value1=reshape(MAdjacency(1:4,1:4),[16,1]);
    
        I1=EntriesLeft1';
        I3=EntriesRight2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(2);
        P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value2=reshape(MAdjacency(1:4,5:8),[16,1]);
    
        I1=EntriesRight2';
        I3=EntriesLeft1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(1);
        P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value3=reshape(MAdjacency(5:8,1:4),[16,1]);
    
        I1=EntriesRight2';
        I3=EntriesRight2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(2);
        P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value4=reshape(MAdjacency(5:8,5:8),[16,1]);
    
        P=[P21;P22;P23;P24];
        P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];


        if length(P)+Indexcounter>length(LatticeRow)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        end
        
        %set the indices
        LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
        LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
        LatticeValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
    
   
    
        Indexcounter=Indexcounter+length(P);
        LatticeCounter=LatticeCounter+1;
    

        %Set the indices for the big matrix
        I1=EntriesRight1';
        I3=EntriesRight1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(1);
        P21=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value1=reshape(MAdjacency(1:4,1:4),[16,1]);
    
        I1=EntriesRight1';
        I3=EntriesLeft2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(1);
        I4=ones(length(I1),1)*currentEdge(2);
        P22=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value2=reshape(MAdjacency(1:4,5:8),[16,1]);
    
        I1=EntriesLeft2';
        I3=EntriesRight1';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(1);
        P23=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value3=reshape(MAdjacency(5:8,1:4),[16,1]);
    
        I1=EntriesLeft2';
        I3=EntriesLeft2';
        I1=repmat(I1,4,1);
        I3=repelem(I3,4,1);
        I2=ones(length(I1),1)*currentEdge(2);
        I4=ones(length(I1),1)*currentEdge(2);
        P24=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        P2Value4=reshape(MAdjacency(5:8,5:8),[16,1]);
    
    
    
        P=[P21;P22;P23;P24];
        P2Value=[P2Value1;P2Value2;P2Value3;P2Value4];

        if length(P)+Indexcounter>length(LatticeRow)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        end
        
        %set the indices
        LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
        LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
        LatticeValue(Indexcounter+1:Indexcounter+length(P))=P2Value;
        
        Indexcounter=Indexcounter+length(P);
        LatticeCounter=LatticeCounter+1;
    end


end

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Edge subdivision was computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%-------------------------------------------------------------------------
%       Block Subdivision
%-------------------------------------------------------------------------

%compute the local used subdivision and adjacency matrices
MAdjacencyTrap3=computeTrapezohedronAdjacencyMatrix(3);
MTrap3=computeTriQuadraticSubdivisionMatrixV2(MAdjacencyTrap3,'PreventInputCheck');

%goes through each block
for i=1:AmountOfBlocks

    %Amount of faces the block lays on
    n=BlockFaces(i);

    %if its 3, use the standard case
    if n == 3
        MAdjacency=MAdjacencyTrap3;
        M=MTrap3;

    %if not, compute the necessary subdivision matrix
    else
        MAdjacency=computeTrapezohedronAdjacencyMatrix(n);
        M=computeTriQuadraticSubdivisionMatrixV2(MAdjacency,'PreventInputCheck');
    end

    %compute indice values
    EntriesRow=zeros(2*n+2,1);
    EntriesCol=zeros(2*n+2,1);
    EntriesRow(1)=1+2*(n)+1;
    EntriesCol(1)=1;

    for j=1:n
        EntriesRow(2+2*(j-1))=1+2*(n)+1+4*(n)+2*(j-1)+1;
        EntriesCol(2+2*(j-1))=1+2*(j-1)+1;

        EntriesRow(2+2*(j-1)+1)=1+2*(n)+1+4*(n)+2*(j-1)+2;
        EntriesCol(2+2*(j-1)+1)=1+2*(j-1)+2;
    end


    EntriesRow(end)=1+2*(n)+1+4*(n)+2*n+1;
    EntriesCol(end)=1+2*(n)+1;


    %compute indices for the big matrix
    I1=EntriesRow;
    I3=EntriesCol; 
    I1=repmat(I1,length(EntriesRow),1);
    I3=repelem(I3,length(EntriesCol),1);
    I2=ones(length(I1),1)*i;
    I4=ones(length(I1),1)*i;

    P2=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P=P2;

    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end
    
    %set the indices
    SBigRow(Indexcounter+1:Indexcounter+length(P))=P;
    SBigCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1);
    SBigValue(Indexcounter+1:Indexcounter+length(P))=reshape(M,[length(MAdjacency)^2,1]);
    
    Indexcounter=Indexcounter+length(P);


    %compute indices for the big matrix
    I1=EntriesRow;
    I3=EntriesRow;
    I1=repmat(I1,length(EntriesRow),1);
    I3=repelem(I3,length(EntriesRow),1);
    I2=ones(length(I1),1)*i;
    I4=ones(length(I1),1)*i;

    P2=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P=P2;

    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end
    
    %set the indices
    LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
    LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
    LatticeValue(Indexcounter+1:Indexcounter+length(P))=reshape(MAdjacency,[length(MAdjacency)^2,1]);

    Indexcounter=Indexcounter+length(P);
    LatticeCounter=LatticeCounter+1;

    %compute indices for the big matrix
    I1=EntriesCol;
    I3=EntriesCol;
    I1=repmat(I1,length(EntriesCol),1);
    I3=repelem(I3,length(EntriesCol),1);
    I2=ones(length(I1),1)*i;
    I4=ones(length(I1),1)*i;

    P2=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
    P=P2;

    if length(P)+Indexcounter>length(LatticeRow)
        [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
    end
    
    %set the indices
    LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
    LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
    LatticeValue(Indexcounter+1:Indexcounter+length(P))=reshape(MAdjacency,[length(MAdjacency)^2,1]);
    
    Indexcounter=Indexcounter+length(P);
    LatticeCounter=LatticeCounter+1;

    %adjacency matrix for the regular part
    MAdjacencyCube=computeTrapezohedronAdjacencyMatrix(3);

    %goes thourh all parts of the block
    for j=1:n

        %set entries
        Entries=zeros(8,1);
        Entries(8)=1+2*(n)+1+4*(j-1)+1;

        EntryMiddle=1+2*(j-1)+1;
        if EntryMiddle == 2
            Entries(2)=1+2*(n);
            Entries(3)=EntryMiddle;
            Entries(4)=EntryMiddle+1;

            Entries(5)=Entries(8)+1;
            Entries(6)=1+2*(n)+1+4*(n)+2*(j-1)+1;
            Entries(7)=1+2*(n)+1+4*(n);
        else
            Entries(2)=EntryMiddle-1;
            Entries(3)=EntryMiddle;
            Entries(4)=EntryMiddle+1;

            Entries(5)=Entries(8)+1;
            Entries(6)=1+2*(n)+1+4*(n)+2*(j-1)+1;
            Entries(7)=Entries(8)-1;
        end

        Entries(1)=1+2*(n)+1;

        %compute indices for the big matrix
        I1=Entries;
        I3=Entries;
        I1=repmat(I1,length(Entries),1);
        I3=repelem(I3,length(Entries),1);
        I2=ones(length(I1),1)*i;
        I4=ones(length(I1),1)*i;
    
        P2=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        
        P=P2;

        if length(P)+Indexcounter>length(LatticeRow)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        end
        
        %set the indices
        LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
        LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
        LatticeValue(Indexcounter+1:Indexcounter+length(P))=reshape(MAdjacencyCube,[length(MAdjacencyCube)^2,1]);
            
        Indexcounter=Indexcounter+length(P);
        LatticeCounter=LatticeCounter+1;


        %set entries
        Entries=zeros(8,1);
        Entries(8)=1+2*(n)+1+4*(j-1)+3;

        EntryMiddle=1+2*(n)+1+4*(n)+2*(j-1)+2;
        if EntryMiddle == 8*n+2
            Entries(2)=EntryMiddle-1;
            Entries(3)=EntryMiddle;
            Entries(4)=1+2*(n)+1+4*(n)+1;

            Entries(5)=Entries(8)+1;
            Entries(6)=1+2*(n);
            Entries(7)=Entries(8)-1;
        else
            Entries(2)=EntryMiddle-1;
            Entries(3)=EntryMiddle;
            Entries(4)=EntryMiddle+1;

            Entries(5)=Entries(8)+1;
            Entries(6)=1+2*(j-1)+2;
            Entries(7)=Entries(8)-1;

        end

        Entries(1)=1+2*(n)+1;


        %compute indices for the big matrix
        I1=Entries;
        I3=Entries;
        I1=repmat(I1,length(Entries),1);
        I3=repelem(I3,length(Entries),1);
        I2=ones(length(I1),1)*i;
        I4=ones(length(I1),1)*i;
    
        P2=computeIndex4D(I1,I2,I3,I4,max(BlockSize),AmountOfBlocks,max(BlockSize));
        
        P=P2;

        if length(P)+Indexcounter>length(LatticeRow)
            [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P);
        end
        
        %set the entries
        LatticeRow(Indexcounter+1:Indexcounter+length(P))=P;
        LatticeCol(Indexcounter+1:Indexcounter+length(P))=ones(length(P),1)*LatticeCounter;
        LatticeValue(Indexcounter+1:Indexcounter+length(P))=reshape(MAdjacencyCube,[length(MAdjacencyCube)^2,1]);
            
        Indexcounter=Indexcounter+length(P);
        LatticeCounter=LatticeCounter+1;
            


    end


end

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Block subdivision was computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%remove empty entries
LatticeCol(LatticeRow==0)=[];
LatticeValue(LatticeRow==0)=[];
LatticeRow(LatticeRow==0)=[];

SBigCol(SBigRow==0)=[];
SBigValue(SBigRow==0)=[];
SBigRow(SBigRow==0)=[];

%Set the Subdivision matrix
SBigSparse=sparse(SBigRow,SBigCol,SBigValue,max(BlockSize)*AmountOfBlocks*max(BlockSize)*AmountOfBlocks,1);

%adjust the size of the subdivision matrix
S=reshape(SBigSparse,[max(BlockSize)*AmountOfBlocks,max(BlockSize)*AmountOfBlocks]);

%removes rows and cols of unnecessary points
Keep=find(max(abs(S),[],2));
S=S(Keep,Keep);

%make a full matrix out of it
S=full(S);

%Set the lattice matrix 
Lattice=sparse(LatticeRow,LatticeCol,LatticeValue,max(BlockSize)*AmountOfBlocks*max(BlockSize)*AmountOfBlocks,AmountOfVolumesBig);

%remove unnecessary entries
KeepLattice=Keep+(Keep'-1)*max(BlockSize)*AmountOfBlocks;
Lattice=Lattice(KeepLattice,:);


if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Subdivision matrix was computed.\n'];
    fprintf(StatusString);
end


%if a visualization is whished we plot the eigenshell of the subdivision
%matrix
if visualization

    %compute the eigenvalues and -vectors
    [V,D]=eigs(S);
    D=diag(D);

    %get the eigenvalues and vectors for 1/2
    D=D([2,3,4]);
    V=V(:,[2,3,4]);
    shouldPrint=true;

    %if they are not the right eigenvalues
    if max(abs(D-0.5))>10^(-8) || max(max(abs(imag(V)))) > 10^(-8)

        %try again with all eigenvalues
        [V,D]=eig(S);
        [Entries,~]=find(abs(D-0.5)<10^(-8));
        D=diag(D);
        D=D(Entries);
        V=V(:,Entries);

        %if they are not the right eigenvalues
        if max(abs(D-0.5))>10^(-8) || max(max(abs(imag(V)))) > 10^(-8)
            B=S-1/2*eye(length(S));
            B=B'*B;
            [V,D]=eig(B);
            D=diag(D);
            a=abs(D)<10^(-10);
            V=V(:,a);
            D=D(a)+1/2;
        end

        %if something went wrong...
        if length(D)<3 || max(abs(D-0.5))>10^(-8)
            warning('Subdominant eigenvalues differ more then 10^(-8) from 1/2. Either there is a problem with the algorithm, or the accuracy is not that high.')
            shouldPrint=false;
        end

    end
    
    %plot the eigenshell
    if shouldPrint

        %orthonormal the Eigenvectors such that the shape is not distorted
        V=orth(V);

        if max(max(abs(imag(V)))) > 10^(-8)
            warning('Eigenvalues have a complex part greater then 10^(-8) due to the eig or eigs function. Visualization might not be correct.')
        end

        V=real(V);

        plotTriQuadraticBSplineLattice(V,Lattice);

        axis equal
        axis off
    
        %adjust the view
        view(View);
    
        %set the title
        title("Eigenshell")
    
        %if the graphic should be printed...
        if PrintImages
            exportgraphics(gca,Prefix+"Eigenshell"+Suffix);
        end
    end
end


end




%Also a local function in computeTriQuadraticSubdivisionMatrixV1Big
function [LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue] = allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P)
% ALLOCATEMORESPACE enlarge the six variables so that they can filled with
% values.
%
% Input: The six variables:             LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue
%        The size of the new values:    P
%
% Output: The six enlarged variables:   LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue
%
% allocateMoreSpace(LatticeRow,LatticeCol,LatticeValue,SBigRow,SBigCol,SBigValue,P)
% enlarge the six variables so that they can filled with values.

%The current size of the variables
Indexcounter=length(LatticeRow);

%save the current values
LatticeRowT=LatticeRow;
LatticeColT=LatticeCol;
LatticeValueT=LatticeValue;

SBigRowT=SBigRow;
SBigColT=SBigCol;
SBigValueT=SBigValue;

%create new larger variables
LatticeRow=zeros(length(LatticeRow)+100000+length(P),1);
LatticeCol=zeros(length(LatticeCol)+100000+length(P),1);
LatticeValue=zeros(length(LatticeValue)+100000+length(P),1);

SBigRow=zeros(length(SBigRow)+100000+length(P),1);
SBigCol=zeros(length(SBigCol)+100000+length(P),1);
SBigValue=zeros(length(SBigValue)+100000+length(P),1);

%set the old values to the new variables
LatticeRow(1:Indexcounter)=LatticeRowT;
LatticeCol(1:Indexcounter)=LatticeColT;
LatticeValue(1:Indexcounter)=LatticeValueT;

SBigRow(1:Indexcounter)=SBigRowT;
SBigCol(1:Indexcounter)=SBigColT;
SBigValue(1:Indexcounter)=SBigValueT;

end









%Also a local function in computeTriQuadraticSubdivisionMatrixV1Big
function [index] = computeIndex4D(i1,i2,i3,i4,m1,m2,m3)
%COMPUTEINDEX4D compute the index for the rows and cols of the big
%tri-quadratic subdivision matrix
%
% Input:    i1 - i4 indices
%           m1 - m3 factors
%
% Output:   index: the corresponding index
%
%computeIndex4D(i1,i2,i3,i4,m1,m2,m3) compute the index for the rows and 
%cols of the big tri-quadratic subdivision matrix

index=i1+m1*(i2-1)+m1*m2*(i3-1)+m1*m2*m3*(i4-1);

end