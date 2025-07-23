function [] = plotTriQuadraticBSplineLattice(Controlpoints,Lattice,varargin)
% PLOTTRIQUADRATICBSPLINELATTICE plots the plotable area of a given
% quadratic B-Spline lattice.
%
% Input:    Controlpoints:  The controlpoints of the lattice
%           Lattice:        The information of geometry. Each column is an
%                           adjacency matrix of one volume of the lattice.
%
% plotTriQuadraticBSplineLattice(Controlpoints,Lattice) plots the plotable 
% area of a given quadratic B-Spline lattice.
%
% plotTriQuadraticBSplineLattice( ,'PointsPerCube',n) plots the B-Spline
% lattice with nxnxn points per cube.
%
%plotTriQuadraticBSplineLattice( ,'ColorTable',CT) plots the B-Spline
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

%The amount of colors
[AmountOfColors,~]=size(ColorTable);

%Adjacency matrix for a cube
CubicAdjacency=[0,1,0,1;...
                1,0,1,0;...
                0,1,0,1;...
                1,0,1,0];

CubicAdjacency=[CubicAdjacency,eye(4);eye(4),CubicAdjacency];

%Graph of the cubic adjacency matrix
CubicGraph=graph(CubicAdjacency);

%computes an adjacency matrix for 2x2x2 cubes. This structure of control
%points can be plotted to a B-Spline cube.
PNew=zeros(27);
PNew([1,4,13,10,2,5,14,11],[1,4,13,10,2,5,14,11])=CubicAdjacency;
CubicAdjacencyBig(:,:,1)=PNew;

PNew=zeros(27);
PNew([7,16,13,4,8,17,14,5],[7,16,13,4,8,17,14,5])=CubicAdjacency;
CubicAdjacencyBig(:,:,2)=PNew;

PNew=zeros(27);
PNew([25,22,13,16,26,23,14,17],[25,22,13,16,26,23,14,17])=CubicAdjacency;
CubicAdjacencyBig(:,:,3)=PNew;

PNew=zeros(27);
PNew([19,10,13,22,20,11,14,23],[19,10,13,22,20,11,14,23])=CubicAdjacency;
CubicAdjacencyBig(:,:,4)=PNew;

PNew=zeros(27);
PNew([3,6,15,12,2,5,14,11],[3,6,15,12,2,5,14,11])=CubicAdjacency;
CubicAdjacencyBig(:,:,5)=PNew;

PNew=zeros(27);
PNew([9,18,15,6,8,17,14,5],[9,18,15,6,8,17,14,5])=CubicAdjacency;
CubicAdjacencyBig(:,:,6)=PNew;

PNew=zeros(27);
PNew([27,24,15,18,26,23,14,17],[27,24,15,18,26,23,14,17])=CubicAdjacency;
CubicAdjacencyBig(:,:,7)=PNew;

PNew=zeros(27);
PNew([21,12,15,24,20,11,14,23],[21,12,15,24,20,11,14,23])=CubicAdjacency;
CubicAdjacencyBig(:,:,8)=PNew;

CubicAdjacencyBig=sum(CubicAdjacencyBig,3);

CubicAdjacencyBig(CubicAdjacencyBig>0)=1;

CubicAdjacencyBigGraph=graph(CubicAdjacencyBig);

%The evaluation points in the domain
EvaluationPoints=0:1/(PointsPerCube-1):1;

%Matrix that converges polynoms to B-Splines
GGG=[1/2,-1,1/2;...
    1/2,1,-1;
    0,0,1/2];

%B-Spline Evaluation points
EvaluationPoints=(EvaluationPoints').^(0:2)*GGG';

%Evaluation matrix
M=computeEvaluationpoints3D(eye(27),EvaluationPoints,EvaluationPoints,EvaluationPoints,3);

%Get the information, which controlpoint is an inner controlpoint. These
%are the candidates for the 2x2x2 structure, which can be plotted
[VS,~,~,~,~,PVertex,IsInnerVertex,~]=RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,1/4);

%computes, which points are active on which volume in the lattice
[AmountOfPoints,~]=size(Controlpoints);
[~,AmountOfVolumes]=size(Lattice);
ActivePoints=zeros(AmountOfPoints,AmountOfVolumes);

for i=1:AmountOfVolumes
    ActivePointsLocal=sum(reshape(Lattice(:,i),[AmountOfPoints,AmountOfPoints]),2);
    ActivePointsLocal(ActivePointsLocal>0)=1;
    ActivePoints(:,i)=ActivePointsLocal;
end


%A counter for the inner vertices (whichs structure arround could be
%possibly plotted)
PVertexCounter=0;

%counter for the color of the cubes
colorCounter=1;


%Variables for the final surf
PlotPointsAGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube);
PlotPointsBGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube);
PlotPointsCGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube);

ColorGlobal=NaN((PointsPerCube+1)*AmountOfVolumes*6,PointsPerCube,3);

NaNCounter=1;


for i=1:AmountOfPoints

    %check several properties. If all are fulfilled, a cube can be plotted

    %Fist condition, the vertex has to be an inner vertex
    if IsInnerVertex(i)

        PVertexCounter=PVertexCounter+1;
        ActiveVolumes=find(ActivePoints(i,:));

        %Second condition: The amount of volumes arround the inner vertex
        %has to be 8
        if length(ActiveVolumes)==8

            %Now check, if the structure is 2x2x2
            Cubic=true;
            for j=1:8
                if Cubic

                    %all volumes around the vertex has to have 8 neighbours
                    LocalP=reshape(Lattice(:,ActiveVolumes(j)),[AmountOfPoints,AmountOfPoints]);
                    ActiveRows=sum(LocalP)>0;
                    if sum(ActiveRows)~=8
                        Cubic=false;
                    else

                        %All volumes arround the vertex has to be cubes.
                        LocalSmallP=LocalP(ActiveRows,ActiveRows);
                        LocalGraph=graph(LocalSmallP);
                        I=isomorphism(LocalGraph,CubicGraph);
                        if isempty(I)
                            Cubic=false;
                        end
                    end
                end
            end

            %if all 8 volumes arround the inner vertex are cubes
            if Cubic

                %check, if the structure is 2x2x2
                LocalPVertex=reshape(PVertex(:,PVertexCounter),[length(VS),length(VS)]);

                [row,~]=find(LocalPVertex);
                row=unique(row);
                LocalPVertex=LocalPVertex(row,row);
                LocalPVertexGraph=graph(LocalPVertex);
                I=isomorphism(LocalPVertexGraph,CubicGraph);

                %if the structure is 2x2x2, the volumes can be plotted
                if ~isempty(I)

                    %compute the structure ov 3x3x3 control points
                    BigAdjacencyMatrix=reshape(sum(Lattice(:,ActiveVolumes),2),[AmountOfPoints,AmountOfPoints]);
                    BigAdjacencyMatrix(BigAdjacencyMatrix>0)=1;
                    [row,~]=find(BigAdjacencyMatrix);
                    row=unique(row);
                    BigAdjacencyMatrix=BigAdjacencyMatrix(row,row);
                    LocalVertices=Controlpoints(row,:);
                    BigAdjacencyMatrixGraph=graph(BigAdjacencyMatrix);
                    I=isomorphism(CubicAdjacencyBigGraph,BigAdjacencyMatrixGraph);
                    LocalVertices=LocalVertices(I,:);
                    
                    %compute the concrete plot points
                    PlotPoints=M*LocalVertices;

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
    end
end


%Plots finally the structure
surf(PlotPointsAGlobal,PlotPointsBGlobal,PlotPointsCGlobal,ColorGlobal)
axis equal
axis off


end
















%Also a local function in plotTriCubicBSplineLattice
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

