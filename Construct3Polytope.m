function [PointCoordinates3D,KiteAll,StatusString,statusCounter] = Construct3Polytope(Input,varargin)
%CONSTRUCT3POLYTOPE computes the 3-polytope to the given graph coded in
% Input. The properties of the polytopes are described in theorem 3.25 of
% the corresponding dissertation.
%
% Input:                adjacency matrix or face matrix or edge matrix of the
%                       polytope
%
% varargin:             Several optional commands (see the manual)
%
% PointCoordinates3D:   The point coordinates of the polytope, of the tangent
%                       points and of the dual polytope.
%
% KiteAll:              The number of the vertices of the kites of the 2D drawing 
%
% StatusString:         The string containing the status information if needed.
%
% statusCounter:        The number of status lines.
%
% For more information see the corresponding manual.



%-------------------------------------------------------------------------
%       Input Parsing
%-------------------------------------------------------------------------


%The amount of additional input parameters  
NumberOfAdditionalInput = nargin - 1;

%Standard values for several options of varargin
status=false;
StatusString="";
statusCounter=1;
visualization=false;
PreventInputCheck=false;
Tolerance=10^(-13);
MaxIterations=100;
StartPointTrials=10;
DrawAfterXKites=1;
PrintImages=false;
Prefix='';
Suffix='';
View=[0,20];

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
    %Specifies whether an Inputcheck should be done 
    elseif ischar(varargin{i}) && strcmp(varargin{i},'PreventInputCheck')
        PreventInputCheck=true;
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-1;
    %Specifies the tolerance of the computation
    elseif ischar(varargin{i}) && strcmp(varargin{i},'Tolerance')
        Tolerance=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;
    %Specifies the maximum amount if iterations of numeric parts
    elseif ischar(varargin{i}) && strcmp(varargin{i},'MaxIterations')
        MaxIterations=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;
        if MaxIterations < 1
            MaxIterations=100;
            warning('MaxIterations has to be a positive number. Set the value to 100.')
        end
    %Specifies the illustration of the quad graph
    elseif ischar(varargin{i}) && strcmp(varargin{i},'DrawAfterXKites')
        DrawAfterXKites=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;
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



%Free line for a better visualization of the status output
if status
    clc 
    StatusString=StatusString+"Compute the subdivision matrix for the given input:\n";
    fprintf(StatusString);
end



%Checks whether the input is a face matrix, an adjacancy matrix or an edge matrix 
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

    
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Edge-List input type was identified.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end

%If the Input is square and every entry is either 0 or 1, then we have an adjacency matrix 
elseif rows==cols && sum(sum((Input==0)+(Input==1)))==rows^2
    AdjacencyMatrix=Input;
    
    %The amount of vertices
    AmountOfVertices=length(AdjacencyMatrix);

    if sum(sum(AdjacencyMatrix==AdjacencyMatrix'))~=rows^2
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
    %The amount of vertices
    AmountOfVertices=max(max(Input));
    
    %Comparing the labeling of the vertices with 1,2,...,AmountOfVertices
    if sum((1:AmountOfVertices)==(unique(Input(Input>0)))')~=AmountOfVertices
        error('Vertices are not named uniformly.');
    end

    
    AdjacencyMatrix=FaceToAdjacencyMatrix(Input);

    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Face matrix input type was identified.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end
end


%-------------------------------------------------------------------------
%       Input Check 
%-------------------------------------------------------------------------

%The user can decide, if the input check sould be done or not. (Advantages
%of preventing the input check is speed)
if ~PreventInputCheck

    %Check for self loops
    if sum(diag(AdjacencyMatrix)==zeros(AmountOfVertices,1))~=AmountOfVertices
        error('Self loops are permited.')
    end
    
    %Check for at least 4 vertices
    if AmountOfVertices<4
        error('A volume needs at least 4 vertices')
    end
    
    %Check for at least 4 edges
    if sum(sum(AdjacencyMatrix))<8
        error('A volume needs at least 4 edges')
    end
    
    
    %Check for 3-vertex-connectivity
    if ~is3Connected(AdjacencyMatrix)
        error('Polytope is not 3-connected.');
    end
    
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Polytope is 3-connected.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end

    %Check for planarity
    if ~PlanarityTest(AdjacencyMatrix)
        error('Polytope is not planar.');
    end
    
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Polytope is planar.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end

else
    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Input check was prevented.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end
end


%-------------------------------------------------------------------------
%       Build up Data Structure
%-------------------------------------------------------------------------

%In the quad graph we have a vertex for every vertex, edge and face of the
%original primal graph. Therefore the labeling in each graph component is
%the following:
%
% 1 : |V| are the name of vertices, which correspond to vertices
% |V|+1 : |V|+|E| are the name of vertices, which correspond to edges
% |V|+|E|+1 : |V|+|E|+|F| are the name of vertices, which correspond to
% faces


PrimalGraph=graph(AdjacencyMatrix);


%Polyeder theorem of Euler
AmountOfFaces=2-PrimalGraph.numnodes+PrimalGraph.numedges;
PrimalEdges=PrimalGraph.Edges.EndNodes;
AmountOfEdges=length(PrimalEdges);

%List of Faces of the Polytope
FaceMatrix=computeFaceMatrix(AdjacencyMatrix,false);


if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Faces were identified.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end


%Building the dual graph and the quad graph
DualGraph=graph([]);
QuadGraph=graph([]);
for i=1:length(PrimalEdges)
    %Find the 2 faces, that have one edge in common.
    [row1,~]=find(FaceMatrix==PrimalEdges(i,1));
    [row2,~]=find(FaceMatrix==PrimalEdges(i,2));
    Connection=row1(ismember(row1,row2));
    %Add a line between these 2 faces which are the nodes in the dual graph
    DualGraph=addedge(DualGraph,Connection(1)+AmountOfVertices+AmountOfEdges,Connection(2)+AmountOfVertices+AmountOfEdges,0);

    %Add 4 lines to the quad graph. 

    % PrimalVertex1 to Edge
    QuadGraph=addedge(QuadGraph,PrimalEdges(i,1),i+AmountOfVertices,0);
    % PrimalVertex2 to Edge
    QuadGraph=addedge(QuadGraph,PrimalEdges(i,2),i+AmountOfVertices,0);
    % DualVertex1 to Edge
    QuadGraph=addedge(QuadGraph,Connection(1)+AmountOfVertices+AmountOfEdges,i+AmountOfVertices,0);
    % DualVertex2 to Edge
    QuadGraph=addedge(QuadGraph,Connection(2)+AmountOfVertices+AmountOfEdges,i+AmountOfVertices,0);

end

AQuadGraph=QuadGraph.adjacency;
Kites=computeFaceMatrix(AQuadGraph,false);

%All Kites are needed for the last step of the polytope drawing
KiteAll=Kites;

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Kitelist was build.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end


%-------------------------------------------------------------------------
%       Compute the Radii of the Quad Graph -> Algorithm 6 of Dissertation
%-------------------------------------------------------------------------


%For the quad graph, one edge of the primal graph has to be deleted.
PrimalDeletedEdge=PrimalEdges(1,:);

%and a corresponding dual edge, which also has to be deleted.
[row1,~]=find(FaceMatrix==PrimalDeletedEdge(1));
[row2,~]=find(FaceMatrix==PrimalDeletedEdge(2));
DualDeletedEdge=row1(ismember(row1,row2))'+AmountOfVertices+AmountOfEdges;

%We now prepare the system of rhos explained in Zie04. Therefore we just
%need the kites, which are inside of the quad graph. In the next step the
%kites outside are deleted. If the deleted kite has a point, which is not
%part of an deleted edge, this boundary point is at the boundary of the
%quad graph. Therefore its ii value is set to pi and not to 2 pi.


%Set the pi value for all points
Pi=2*pi*ones(AmountOfVertices+AmountOfEdges+AmountOfFaces,1);

%Goes through all Kites and delete the outer (not inside the quad) kites.
[a,~]=size(Kites);
i=1;
while i<=a
    if ismember(min(Kites(i,:)),PrimalDeletedEdge) && ismember(max(Kites(i,:)),DualDeletedEdge)
        Kites(i,:)=[];
        a=a-1;
    elseif ismember(min(Kites(i,:)),PrimalDeletedEdge)

        %If one point is on the boundary, Pi value is reduced
        Pi(max(Kites(i,:)))=pi;
        Kites(i,:)=[];
        a=a-1;
    elseif ismember(max(Kites(i,:)),DualDeletedEdge)

        %If one point is on the boundary, Pi value is reduced
        Pi(min(Kites(i,:)))=pi;
        Kites(i,:)=[];
        a=a-1;
    else
        i=i+1;
    end
end

%Kombination of primal and dual points via kites 
FCombinations=[min(Kites,[],2),max(Kites,[],2)];

%RhoIdentification is the mapping of variables of the system to the points 
[RhoIdentifications,~,RhoCombinations]=unique(FCombinations);

%The rho combinations in terms of uniform named variables. 
RhoCombinations=reshape(RhoCombinations,a,2);

%Pi reduced to the variables of the system
PiBig=Pi;
Pi=Pi(RhoIdentifications);

%-------------------------------------------------------------------------
%Try to solve the system with Newtons method

%Create starting point
EvaluationPoint=zeros(length(Pi),1);

%Evaluate gradient
Gradient=EvaluateSRho(EvaluationPoint,RhoCombinations,Pi,true);

%Iteration counter
Iterations=0;

%Do newtons method as long as the solution is found (better then the given
%tolerance) or the maximum Iteration is reached.
while max(abs(Gradient)) >Tolerance && Iterations < MaxIterations

    %Evaluate Hesse Matrix
    Hesse=EvaluateSRhoDerivative(EvaluationPoint,RhoCombinations,true);

    %Compute Direction
    Direction=Hesse\(Gradient);

    %Set new evaluation point
    EvaluationPoint=EvaluationPoint-Direction;

    %Evaluate new gradient
    Gradient=EvaluateSRho(EvaluationPoint,RhoCombinations,Pi,true);

    %Increase Iterations
    Iterations=Iterations+1;


    if status
        clc 
        fprintf(StatusString);
        disp(['(', num2str(statusCounter) , ') Try to compute radii by Netwons Method'])
        disp(['    Iteration: ', num2str(Iterations),' / ',num2str(MaxIterations)  ])
        disp(['    Error: ', num2str(max(abs(Gradient))),' / ',num2str(Tolerance) ])
    end
end

%If Newtons method was sucessfull, set rho as the gradient
if max(abs(Gradient)) <=Tolerance && sum(isnan(Gradient)) == 0
    rho=EvaluationPoint;
    SolutionFound=true;
%If not use the standard solver of matlab.
else

    if status
        clc 
        StatusString=StatusString+['(', num2str(statusCounter) , ') Computing the radii with Newtons method was not sucessful. Try to use fsolve.\n'];
        fprintf(StatusString);
        statusCounter=statusCounter+1;
    end

    %Set the EvaluationPoint to zero
    EvaluationPointNewton=EvaluationPoint;
    EvaluationPoint=zeros(length(Pi),1);

    %Initiate the gradient function
    f=@(x) EvaluateSRho(x,RhoCombinations,Pi,true);
    
    %Set the options for the standard solver
    options = optimoptions('fsolve','FunctionTolerance',Tolerance,'StepTolerance',Tolerance,'OptimalityTolerance',Tolerance,'Display','off');
    
    %Set the solutionFound value and the try counter
    SolutionFound=false;
    TryCounter=1;
    
    %As long as we do not find a solution or as the maximum of trials was
    %not reached
    while ~SolutionFound && TryCounter < StartPointTrials
        
        if status
            clc 
            fprintf(StatusString);
            disp(['    StartPointTrials: ', num2str(TryCounter),' / ',num2str(StartPointTrials)  ])
        end

        %Use matlab standard solver to solve the system (wanrnings needed
        %to be off as the solver will produce a warning as the system is
        %non quadratic.
        warning off;
        [rho,~,exitflag]=fsolve(f,EvaluationPoint,options);
        warning on;
        
        %If the solver did not find a solution
        if exitflag < 0 && max(abs(EvaluateSRho(rho,RhoCombinations,Pi,true)))> Tolerance

            %Try another starting point
            TryCounter=TryCounter+1;

            %As second try use the gradient computed in the previous method
            if TryCounter == 2
                EvaluationPoint=EvaluationPointNewton;
            else
                EvaluationPoint=rand(length(Pi),1);
            end

        %If the algorithm thinks she found a solution but the solution is not valid    
        elseif max(abs(EvaluateSRho(rho,RhoCombinations,Pi,true))) > Tolerance
            error('fsolve found a solution, but this solution was worse than the given tolerance.')

        %Otherwise the solution was valid
        else
            if sum(isnan(rho))==0
                SolutionFound=true;
            else
                TryCounter=TryCounter+1;
                EvaluationPoint=rand(length(Pi),1);
            end
        end
    end

end

if ~SolutionFound
    error('It was not possible to find a list of radii within the given tolerance.')
end

%Set the radii
r=exp(rho);

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Radii were computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end

%-----------------------------------------------------------------------------
%       Compute the Coordinates of the Quad Graph Algorithm 7 of Dissertation
%-----------------------------------------------------------------------------

%As we have now all length of all lines we compute now the coordinates of
%the points of the quad graph. The idea is that we draw a kite after
%another by a BFS search. If the relative error of the points went to big, 
% we readjust all points.

%Initiate the variables for the BFS search

%The radii of the numbered points (just vertex and face points have a given
%radius).
RPoints=zeros(AmountOfVertices+AmountOfEdges+AmountOfFaces,1);
RPoints(RhoIdentifications)=r;

%The number of kites which has to be drawn
[AmountOfKites,~]=size(Kites);

%A list containing the information if a kite was drawn or not
KiteSeen=zeros(AmountOfKites,1);


%A list containing the information if a point was drawn or not
PointSeen=zeros(AmountOfVertices+AmountOfEdges+AmountOfFaces,1);

%The concrete koordinates of the drawn points 
PointCoordinates=zeros(AmountOfVertices+AmountOfEdges+AmountOfFaces,2);

%A list of all edges
Edges=QuadGraph.Edges.EndNodes;

%The orientation each edge was drawn (for the left, right information on
%which site the new kite should be.
EdgeOrientation=zeros(length(Edges),1);

%A priority list of the edges (decide, which kite should be drawn next).
EdgeOrder=zeros(length(Edges),1)+1000+AmountOfKites;

%Lists of the edges depending of the kites
Edge1=Kites(:,[1,2]);
Edge2=Kites(:,[2,3]);
Edge3=Kites(:,[4,3]);
Edge4=Kites(:,[1,4]);

%Link number of the EdgeKiteLists to the general EdgeList
[~,Edge1Number]=ismember(Edge1,Edges,'rows');
[~,Edge2Number]=ismember(Edge2,Edges,'rows');
[~,Edge3Number]=ismember(Edge3,Edges,'rows');
[~,Edge4Number]=ismember(Edge4,Edges,'rows');


%Searches for an edge kite (a kite for which the primal and dual point are
%both at the boundary. This will be the starting kite of the algorithm.
found=false;
counter=1;
while ~found
    if PiBig(Kites(counter,1))==pi && PiBig(Kites(counter,3))==pi
        found=true;
        CornerKite=counter;
    end
    counter=counter+1;
end

%Set the values for the first kite


%Check which of the edge points lay on the corner
UniqueVertex=find(Kites==Kites(CornerKite,2));

if isscalar(UniqueVertex)
    %The corner point
    PointCoordinates(Kites(CornerKite,2),:)=[0,0];
    
    %The first primal point
    PointCoordinates(Kites(CornerKite,1),:)=[0,-RPoints(Kites(CornerKite,1))];
    
    %The first edge
    CurrentEdge=[min(Kites(CornerKite,1:2),max(Kites(CornerKite,1:2)))];

else
    %The corner point
    PointCoordinates(Kites(CornerKite,4),:)=[0,0];

    %The first primal point
    PointCoordinates(Kites(CornerKite,1),:)=[0,-RPoints(Kites(CornerKite,1))];
    
    %The first edge
    CurrentEdge=[min(Kites(CornerKite,[1,4]),max(Kites(CornerKite,[1,4])))];
end

[~,Position]=ismember(CurrentEdge,Edges,'rows');
EdgeOrientation(Position)=1;
EdgeOrder(Position)=1;

%Compute the restrictions for the kites/triangles
    
%Every primal edge has to have the defined length of the radius of the
%primal point
PointCombinationPrimalEdge=unique([Kites(:,[1,2]);Kites(:,[1,4])],'rows');
rrPrimalEdge=RPoints(PointCombinationPrimalEdge(:,1));

%Every dual edge has to have the defined length of the radius of the
%dual point
PointCombinationDualEdge=unique([Kites(:,[3,2]);Kites(:,[3,4])],'rows');
rrDualEdge=RPoints(PointCombinationDualEdge(:,1));

%Every diagonal has to respect the right angle of the primal dual
%combination
PointCombinationPrimalDual=[Kites(:,[1,3])];
rrPrimalDual=(sqrt(RPoints(Kites(:,1)).^2+RPoints(Kites(:,3)).^2));

%Two triangles in one kite do not have to overlap
PointCombinationEdgeEdge=Kites(:,[2,4]);
rrEdgeEdge=2*RPoints(Kites(:,1)).*RPoints(Kites(:,3))./(sqrt(RPoints(Kites(:,1)).^2+RPoints(Kites(:,3)).^2));


%Two kites does not have to overlap. Thus we check if every primal point
%and every dual point sharing a common edge point have radii+radii distance
EdgePoints=unique(Kites(:,[2,4]));

%A list of all possible combinations
PointCombinationPrimalPrimal=zeros(length(EdgePoints),2);
PointCombinationDualDual=zeros(length(EdgePoints),2);

PrimalCounter=1;
DualCounter=1;

%Goes through all edge point
for i=1:length(EdgePoints)
    [row,~]=find(Kites==EdgePoints(i));
    
    %Check, if an edge point has two primal points
    if min(unique(Kites(row,1)))~=max(unique(Kites(row,1)))

        %If so, add them to the list
        PointCombinationPrimalPrimal(PrimalCounter,:)=[min(unique(Kites(row,1))),max(unique(Kites(row,1)))];
        PrimalCounter=PrimalCounter+1;
    end
    
    %Check, if an edge point has two dual points
    if min(unique(Kites(row,3)))~=max(unique(Kites(row,3)))
        
        %If so, add them to the list
        PointCombinationDualDual(DualCounter,:)=[min(unique(Kites(row,3))),max(unique(Kites(row,3)))];
        DualCounter=DualCounter+1;
    end

end

%Reduce the combination to the real entries
PointCombinationPrimalPrimal=PointCombinationPrimalPrimal(1:PrimalCounter-1,:);

%Unique them to avoid doubles
PointCombinationPrimalPrimal=unique(PointCombinationPrimalPrimal,'rows');

%Compute the right site (radius+radius)
rrPrimalPrimal=RPoints(PointCombinationPrimalPrimal(:,1))+RPoints(PointCombinationPrimalPrimal(:,2));

%Reduce the combination to the real entries
PointCombinationDualDual=PointCombinationDualDual(1:DualCounter-1,:);

%Unique them to avoid doubles
PointCombinationDualDual=unique(PointCombinationDualDual,'rows');

%Compute the right site (radius+radius)
rrDualDual=RPoints(PointCombinationDualDual(:,1))+RPoints(PointCombinationDualDual(:,2));

%Put all restrictions and their right site together
PointCombinationGes=[PointCombinationPrimalEdge;PointCombinationDualEdge;PointCombinationPrimalDual;PointCombinationEdgeEdge;PointCombinationPrimalPrimal;PointCombinationDualDual];
rrGes=[rrPrimalEdge;rrDualEdge;rrPrimalDual;rrEdgeEdge;rrPrimalPrimal;rrDualDual];

%The amount of restrictions
AmountOfRestrictions=length(PointCombinationGes);

%Tells, if one point in a restriction was already seen
SeenInPointCombination=zeros(AmountOfRestrictions,2);

%Tells, if the whole restriction was already seen
EnabledRestriction=false(AmountOfRestrictions,1);

%Resets the current figure if necessary
if visualization && DrawAfterXKites>0
    figure
end




%Starting the BFS algorithm. Go trough all kites by a neighbor relation and
%draw each kite step by step. If the accuracy is too bad, we readjust the
%points by Newtons method
while sum(KiteSeen)<AmountOfKites

    %Status update
    if status
        clc 
        fprintf(StatusString);
        disp(['(', num2str(statusCounter) , ') Drawing the Kites'])
        disp(['    Kite drawn: ', num2str(sum(KiteSeen)),' / ',num2str(AmountOfKites)  ])
    end


    %Searches the next kite, the kite with an edge of the lowest priority

    %Initiate Start values
    CurrentKite=0;
    CurrentKitePriority=1000+AmountOfKites;

    %Goes through all kites
    for i = 1:AmountOfKites

         %If the Kites was not drawn...
         if KiteSeen(i)==0

           LocalKitePriority=min([EdgeOrder(Edge1Number(i)),EdgeOrder(Edge2Number(i)),EdgeOrder(Edge3Number(i)),EdgeOrder(Edge4Number(i))]);

           %...and the priority of those edges is lower then the current
           %found...
           if LocalKitePriority < CurrentKitePriority

                %...then set the kite as the current kite
                CurrentKite=i;
                CurrentKitePriority=LocalKitePriority;
           end
        end
    end
    
    %Set the point identification numbers of the current kite
    PointNumberPrimal=Kites(CurrentKite,1);
    PointNumberE1=Kites(CurrentKite,2);
    PointNumberDual=Kites(CurrentKite,3);
    PointNumberE2=Kites(CurrentKite,4);

    %Set the current minimal edge order (the edge of the kite with the minimal
    %order)
    MinEdgeOrder=min([EdgeOrder(Edge1Number(CurrentKite)),EdgeOrder(Edge2Number(CurrentKite)),EdgeOrder(Edge3Number(CurrentKite)),EdgeOrder(Edge4Number(CurrentKite))]);


    %Drawing of the kite: The drawing depends on the edge which was visited.
    %Therefore we need 4 different cases (one for each edge)
    if EdgeOrder(Edge1Number(CurrentKite)) == MinEdgeOrder

        %-------------------------------------------------------------------------
        % 1. edge

        %Get the coordinates of the points which already be drawn (as
        %endpoints of the edge, through which the kite was visited
        PointKoordinatesPrimal=PointCoordinates(PointNumberPrimal,:);
        PointKoordinatesE1=PointCoordinates(PointNumberE1,:);

        %Compute the direction of the visiting edge
        Direction=-PointKoordinatesE1+PointKoordinatesPrimal;
        Direction=Direction./norm(Direction);
        Direction=Direction';

        %Depending on the orientation of the edge, the direction is rotated
        %clockwise or counterclockwise with 90 degree. The orientation of
        %the other edges are set for the other neighboring kites
        if EdgeOrientation(Edge1Number(CurrentKite))==1

            EdgeOrientation(Edge2Number(CurrentKite))=1;
            EdgeOrientation(Edge4Number(CurrentKite))=1;
            EdgeOrientation(Edge3Number(CurrentKite))=-1;

            Direction=[0,-1;1,0]*Direction;
        else
            EdgeOrientation(Edge2Number(CurrentKite))=-1;
            EdgeOrientation(Edge4Number(CurrentKite))=-1;
            EdgeOrientation(Edge3Number(CurrentKite))=1;

            Direction=[0,1;-1,0]*Direction;
        end
        Direction=Direction';

        %Compute the first new point as the edgepoint plus the computed
        %direction (angle=90 degree) times the given radius.
        PointKoordinatesDual=PointKoordinatesE1+Direction*RPoints(PointNumberDual);

        %Compute the second new point by computing the rectangular
        %projection of e_1 on PD and adding twiche the distance to e_1
        R=[0,-1;1,0];
        V1=(-PointKoordinatesPrimal+PointKoordinatesDual)';
        V2=-R*(-PointKoordinatesPrimal+PointKoordinatesDual)';
        RS=(PointKoordinatesE1-PointKoordinatesPrimal)';
        Parameter=[V1,V2]\RS;
        Intersection=PointKoordinatesPrimal+(-PointKoordinatesPrimal+PointKoordinatesDual)*Parameter(1);
        PointKoordinatesE2=PointKoordinatesE1+2*(-PointKoordinatesE1+Intersection);


    elseif EdgeOrder(Edge2Number(CurrentKite)) == MinEdgeOrder

        %-------------------------------------------------------------------------
        % 2. Edge

        %Get the coordinates of the points which already be drawn (as
        %endpoints of the edge, through which the kite was visited
        PointKoordinatesDual=PointCoordinates(PointNumberDual,:);
        PointKoordinatesE1=PointCoordinates(PointNumberE1,:);
        
        %Compute the direction of the visiting edge
        Direction=-PointKoordinatesE1+PointKoordinatesDual;
        Direction=Direction./norm(Direction);
        Direction=Direction';

        %Depending on the orientation of the edge, the direction is rotated
        %clockwise or counterclockwise with 90 degree. The orientation of
        %the other edges are set for the other neighboring kites
        if EdgeOrientation(Edge2Number(CurrentKite))==1

            EdgeOrientation(Edge1Number(CurrentKite))=1;
            EdgeOrientation(Edge3Number(CurrentKite))=1;
            EdgeOrientation(Edge4Number(CurrentKite))=-1;

            Direction=[0,-1;1,0]*Direction;
        else

            EdgeOrientation(Edge1Number(CurrentKite))=-1;
            EdgeOrientation(Edge3Number(CurrentKite))=-1;
            EdgeOrientation(Edge4Number(CurrentKite))=1;

            Direction=[0,1;-1,0]*Direction;
        end
        Direction=Direction';

        %Compute the first new point as the edgepoint plus the computed
        %direction (angle=90 degree) times the given radius.
        PointKoordinatesPrimal=PointKoordinatesE1+Direction*RPoints(PointNumberPrimal);

        %Compute the second new point by computing the rectangular
        %projection of e_1 on PD and adding twiche the distance to e_1
        R=[0,-1;1,0];
        V1=(-PointKoordinatesPrimal+PointKoordinatesDual)';
        V2=-R*(-PointKoordinatesPrimal+PointKoordinatesDual)';
        RS=(PointKoordinatesE1-PointKoordinatesPrimal)';
        Parameter=[V1,V2]\RS;
        Intersection=PointKoordinatesPrimal+(-PointKoordinatesPrimal+PointKoordinatesDual)*Parameter(1);
        PointKoordinatesE2=PointKoordinatesE1+2*(-PointKoordinatesE1+Intersection);


    elseif EdgeOrder(Edge3Number(CurrentKite)) == MinEdgeOrder

        %-------------------------------------------------------------------------
        % 3. Edge

        %Get the coordinates of the points which already be drawn (as
        %endpoints of the edge, through which the kite was visited
        PointKoordinatesDual=PointCoordinates(PointNumberDual,:);
        PointKoordinatesE2=PointCoordinates(PointNumberE2,:);

        %Compute the direction of the visiting edge
        Direction=-PointKoordinatesE2+PointKoordinatesDual;
        Direction=Direction./norm(Direction);
        Direction=Direction';

        %Depending on the orientation of the edge, the direction is rotated
        %clockwise or counterclockwise with 90 degree. The orientation of
        %the other edges are set for the other neighboring kites
        if EdgeOrientation(Edge3Number(CurrentKite))==1

            EdgeOrientation(Edge2Number(CurrentKite))=1;
            EdgeOrientation(Edge4Number(CurrentKite))=1;
            EdgeOrientation(Edge1Number(CurrentKite))=-1;


            Direction=[0,-1;1,0]*Direction;
        else
            EdgeOrientation(Edge2Number(CurrentKite))=-1;
            EdgeOrientation(Edge4Number(CurrentKite))=-1;
            EdgeOrientation(Edge1Number(CurrentKite))=1;

            Direction=[0,1;-1,0]*Direction;
        end
        Direction=Direction';

        %Compute the first new point as the edgepoint plus the computed
        %direction (angle=90 degree) times the given radius.
        PointKoordinatesPrimal=PointKoordinatesE2+Direction*RPoints(PointNumberPrimal);

        %Compute the second new point by computing the rectangular
        %projection of e_1 on PD and adding twice the distance to e_1
        R=[0,-1;1,0];
        V1=(-PointKoordinatesPrimal+PointKoordinatesDual)';
        V2=-R*(-PointKoordinatesPrimal+PointKoordinatesDual)';
        RS=(PointKoordinatesE2-PointKoordinatesPrimal)';
        Parameter=[V1,V2]\RS;
        Intersection=PointKoordinatesPrimal+(-PointKoordinatesPrimal+PointKoordinatesDual)*Parameter(1);
        PointKoordinatesE1=PointKoordinatesE2+2*(-PointKoordinatesE2+Intersection);


    else
        %-------------------------------------------------------------------------
        % 4. Edge

        %Get the coordinates of the points which already be drawn (as
        %endpoints of the edge, through which the kite was visited
        PointKoordinatesPrimal=PointCoordinates(PointNumberPrimal,:);
        PointKoordinatesE2=PointCoordinates(PointNumberE2,:);
        
        %Compute the direction of the visiting edge
        Direction=-PointKoordinatesE2+PointKoordinatesPrimal;
        Direction=Direction./norm(Direction);
        Direction=Direction';

        %Depending on the orientation of the edge, the direction is rotated
        %clockwise or counterclockwise with 90 degree. The orientation of
        %the other edges are set for the other neighboring kites
        if EdgeOrientation(Edge4Number(CurrentKite))==1

            EdgeOrientation(Edge1Number(CurrentKite))=1;
            EdgeOrientation(Edge3Number(CurrentKite))=1;
            EdgeOrientation(Edge2Number(CurrentKite))=-1;

            Direction=[0,-1;1,0]*Direction;
        else

            EdgeOrientation(Edge1Number(CurrentKite))=-1;
            EdgeOrientation(Edge3Number(CurrentKite))=-1;
            EdgeOrientation(Edge2Number(CurrentKite))=1;

            Direction=[0,1;-1,0]*Direction;
        end
        Direction=Direction';

        %Compute the first new point as the edgepoint plus the computed
        %direction (angle=90 degree) times the given radius.
        PointKoordinatesDual=PointKoordinatesE2+Direction*RPoints(PointNumberDual);

        %Compute the second new point by computing the rectangular
        %projection of e_1 on PD and adding twiche the distance to e_1
        R=[0,-1;1,0];
        V1=(-PointKoordinatesPrimal+PointKoordinatesDual)';
        V2=-R*(-PointKoordinatesPrimal+PointKoordinatesDual)';
        RS=(PointKoordinatesE2-PointKoordinatesPrimal)';
        Parameter=[V1,V2]\RS;
        Intersection=PointKoordinatesPrimal+(-PointKoordinatesPrimal+PointKoordinatesDual)*Parameter(1);
        PointKoordinatesE1=PointKoordinatesE2+2*(-PointKoordinatesE2+Intersection);
    end

    %Set the current kite as seen
    KiteSeen(CurrentKite)=1;

    %Set the edge priority of the other edges to the current priority +1
    EdgeOrder([Edge1Number(CurrentKite),Edge2Number(CurrentKite),Edge3Number(CurrentKite),Edge4Number(CurrentKite)])=CurrentKitePriority+1;
    
    %Check for every drawn point if it was already seen. If so, compute the
    %arithmetic mean of both of them and write the coordinates of them
    if PointSeen(Kites(CurrentKite,1))==0
        PointCoordinates(Kites(CurrentKite,1),:)=PointKoordinatesPrimal;
        PointSeen(Kites(CurrentKite,1))=1;
        [row,col]=find(PointCombinationGes==Kites(CurrentKite,1));
        SeenInPointCombination(row+(col-1)*AmountOfRestrictions)=1;
    else
        PointCoordinates(Kites(CurrentKite,1),:)=(PointCoordinates(Kites(CurrentKite,1),:)+PointKoordinatesPrimal)/2;
    end

    if PointSeen(Kites(CurrentKite,2))==0
        PointCoordinates(Kites(CurrentKite,2),:)=PointKoordinatesE1;
        PointSeen(Kites(CurrentKite,2))=1;
        [row,col]=find(PointCombinationGes==Kites(CurrentKite,2));
        SeenInPointCombination(row+(col-1)*AmountOfRestrictions)=1;
    else
        PointCoordinates(Kites(CurrentKite,2),:)=(PointCoordinates(Kites(CurrentKite,2),:)+PointKoordinatesE1)/2;
    end

    if PointSeen(Kites(CurrentKite,3))==0
        PointCoordinates(Kites(CurrentKite,3),:)=PointKoordinatesDual;
        PointSeen(Kites(CurrentKite,3))=1;
        [row,col]=find(PointCombinationGes==Kites(CurrentKite,3));
        SeenInPointCombination(row+(col-1)*AmountOfRestrictions)=1;
    else
        PointCoordinates(Kites(CurrentKite,3),:)=(PointCoordinates(Kites(CurrentKite,3),:)+PointKoordinatesDual)/2;
    end

    if PointSeen(Kites(CurrentKite,4))==0
        PointCoordinates(Kites(CurrentKite,4),:)=PointKoordinatesE2;
        PointSeen(Kites(CurrentKite,4))=1;
        [row,col]=find(PointCombinationGes==Kites(CurrentKite,4));
        SeenInPointCombination(row+(col-1)*AmountOfRestrictions)=1;
    else
        PointCoordinates(Kites(CurrentKite,4),:)=(PointCoordinates(Kites(CurrentKite,4),:)+PointKoordinatesE2)/2;
    end
    
    %Check which restrictions are active and sets it
    EnabledRestriction(sum(SeenInPointCombination,2)==2)=1;

    %Get just the restrictions, which are active
    PointKombinationLocal=PointCombinationGes(EnabledRestriction,:);
    rrLocal=rrGes(EnabledRestriction);
    
    %Compute the error of the restrictions
    Gradient=EvaluateDistance(PointCoordinates(:,1),PointCoordinates(:,2),PointKombinationLocal(:,:),rrLocal(:));

    %Save the original point coordinates
    PointKoordinatesOld=PointCoordinates;

    %Set the iterations of Newtons Method (to avoid infinity loop
    Iterations=1;

    %As long as at least one restriction is hurt
    while max(abs(Gradient)) > Tolerance*max(rrGes) && Iterations < MaxIterations
        
        %Status update
        if status
            clc 
            fprintf(StatusString);
            disp(['(', num2str(statusCounter) , ') Drawing the Kites'])
            disp(['    Kite drawn: ', num2str(sum(KiteSeen)),' / ',num2str(AmountOfKites)  ])
            disp('    At least one restriction was hurt Rearange points with newtons method.  ')
            disp(['    Iteration: ', num2str(Iterations),' / ',num2str(MaxIterations)  ])
            disp(['    Error: ', num2str(max(abs(Gradient))),' / ',num2str(Tolerance*max(rrGes)) ])
        end

        %Evaluate Hesse-Matrix
        Hesse=EvaluateDistanceDerivative(PointCoordinates(:,1),PointCoordinates(:,2),PointKombinationLocal(:,:));
    
        %Get just the columns of the active points (at first all x
        %coordinates, then all y coordinates)
        HesseLocal=Hesse(:,logical([PointSeen;PointSeen]));

        %Compute direction (as the system might be linear dependend because
        %of redundand restrictions, warnings are off
        warning off
        Direction=HesseLocal\(Gradient);
        warning on
        
        %Set new evaluation point
        PointCoordinates(logical(PointSeen),:)=reshape(reshape(PointCoordinates(logical(PointSeen),:),[2*sum(PointSeen),1])-Direction,[sum(PointSeen),2]);
    
        %Evaluate new gradient
        Gradient=EvaluateDistance(PointCoordinates(:,1),PointCoordinates(:,2),PointKombinationLocal(:,:),rrLocal(:));

        Iterations=Iterations+1;
    end

    %If after the upper Newton's method one restriction is still hurt...
    if max(abs(Gradient)) > Tolerance*max(rrGes) || sum(isnan(Gradient))>0

        if status
            clc 
            fprintf(StatusString);
            disp(['(', num2str(statusCounter) , ') Drawing the Kites'])
            disp(['    Kite drawn: ', num2str(sum(KiteSeen)),' / ',num2str(AmountOfKites)  ])
            disp('    At least one restriction was hurt Rearange points with newtons method.' )
            disp('    Newtons method was not sucessful. Try fsolve.')
            
        end

        %...we use fsolve with the original points
        PointCoordinates=PointKoordinatesOld;
        f=@(x) EvaluateDistance(x(:,1),x(:,2),PointKombinationLocal(:,:),rrLocal(:));

        %As we do not have a quadratic system warnings are off
        warning off
        options = optimoptions('fsolve','FunctionTolerance',Tolerance*max(rrGes),'StepTolerance',Tolerance*max(rrGes),'OptimalityTolerance',Tolerance*max(rrGes),'Display','off');
        [PointCoordinates,~,]=fsolve(f,PointCoordinates,options); %(PointSeen,:)
        warning on

        %If a restiction is still hurt we throw an error
        Gradient=EvaluateDistance(PointCoordinates(:,1),PointCoordinates(:,2),PointKombinationLocal(:,:),rrLocal(:));
        if max(abs(Gradient)) > Tolerance*max(rrGes) || sum(isnan(Gradient))>0
            error('Could not draw the quad graph within the given Tolerance.')
        end

    end
    
    %plots the construction steps of the quad graph
    if visualization
        hold off

        %the first step is just a line and therefore different
        if sum(KiteSeen)==1
           
            %goes through all kites...
            for i = 1:AmountOfKites
                    
                %and searches for the first seen one
                if KiteSeen(i)==1

                    %The first edge is an edge between an edge point and a
                    %primal point. In the beginning it is the only edge,
                    %with a priority, therefore it is plotted red
                    plot(PointCoordinates(Kites(i,[1,2]),1),PointCoordinates(Kites(i,[1,2]),2),'Color',[1,0,0],'LineWidth',3);
                    hold on
                    axis equal

                    %plots the priority number of the red edge
                    V=[0,-1;1,0]*[diff(PointCoordinates(Kites(i,[1,2]),1));diff(PointCoordinates(Kites(i,[1,2]),2))];
                    V=-V/norm(V);
                    text(sum(PointCoordinates(Kites(i,[1,2]),1))/2+V(1)*0.1,sum(PointCoordinates(Kites(i,[1,2]),2))/2+V(2)*0.1,int2str( 1))
                  
                    %plots the primal and edge point
                    plot(PointCoordinates(Kites(i,1),1),PointCoordinates(Kites(i,1),2),'.','Color',[0.5,0,0],'Markersize',40);
                    plot(PointCoordinates(Kites(i,2),1),PointCoordinates(Kites(i,2),2),'.','Color',[0,0.5,0],'Markersize',40);
    
                end
            end

            %set the title
            title("Quad Graph")
    
            %if the image should printed it is done here
            if PrintImages
                exportgraphics(gca,Prefix+"QuadGraphBuild"+int2str((sum(KiteSeen)-1))+Suffix);
            end

        end

        %plots the current step of the drawing process (only if it is in
        %the right ireration described by DrawAfterXKites
        if mod(sum(KiteSeen),DrawAfterXKites)==0

            hold off

            %variable, that activates graphical commands after the first
            % point is plotted
            FirstPoint=true;
            for i = 1:AmountOfKites
                
                %if the kite is already drawn
                if KiteSeen(i)==1
                    

                   %draws one kite after another. Each line is checked, if
                   %all adjacent lines are alredy drawn (black line) or not
                   %(red line)

                   %Checks if red or black
                   [rows1,~]=find(Kites==Kites(i,1));
                   [rows2,~]=find(Kites==Kites(i,2));
                   AK=rows1(ismember(rows1,rows2));
                   if sum(KiteSeen(AK))~=length(AK)

                        %if red, the priority number is calculated and
                        %plotted
                        [~,EN]=ismember([min(Kites(i,[1,2])),max(Kites(i,[1,2]))],Edges,'rows');
                        plot(PointCoordinates(Kites(i,[1,2]),1),PointCoordinates(Kites(i,[1,2]),2),'Color',[1,0,0],'LineWidth',3);
                        V=[0,-1;1,0]*[diff(PointCoordinates(Kites(i,[1,2]),1));diff(PointCoordinates(Kites(i,[1,2]),2))];
                        V=V/norm(V);
                        text(sum(PointCoordinates(Kites(i,[1,2]),1))/2+V(1)*0.1,sum(PointCoordinates(Kites(i,[1,2]),2))/2+V(2)*0.1,int2str( EdgeOrder(EN)))
                   else
                        plot(PointCoordinates(Kites(i,[1,2]),1),PointCoordinates(Kites(i,[1,2]),2),'Color','k','LineWidth',3);
                   end
                
                   if FirstPoint
                        hold on
                        axis equal
                        FirstPoint=false;
                   end
                   
                   
                   %Checks if red or black
                   [rows1,~]=find(Kites==Kites(i,2));
                   [rows2,~]=find(Kites==Kites(i,3));
                   AK=rows1(ismember(rows1,rows2));
                   if sum(KiteSeen(AK))~=length(AK)

                       %if red, the priority number is calculated and
                       %plotted
                       [~,EN]=ismember([min(Kites(i,[2,3])),max(Kites(i,[2,3]))],Edges,'rows');
                       plot(PointCoordinates(Kites(i,[2,3]),1),PointCoordinates(Kites(i,[2,3]),2),'Color',[1,0,0],'LineWidth',3);
                       V=[0,-1;1,0]*[diff(PointCoordinates(Kites(i,[2,3]),1));diff(PointCoordinates(Kites(i,[2,3]),2))];
                       V=V/norm(V);
                       text(sum(PointCoordinates(Kites(i,[2,3]),1))/2+V(1)*0.1,sum(PointCoordinates(Kites(i,[2,3]),2))/2+V(2)*0.1,int2str( EdgeOrder(EN)))
                   else
                       plot(PointCoordinates(Kites(i,[2,3]),1),PointCoordinates(Kites(i,[2,3]),2),'Color','k','LineWidth',3);
                   end

                   
                   %Checks if red or black
                   [rows1,~]=find(Kites==Kites(i,3));
                   [rows2,~]=find(Kites==Kites(i,4));
                   AK=rows1(ismember(rows1,rows2));
                   if sum(KiteSeen(AK))~=length(AK)

                       %if red, the priority number is calculated and
                       %plotted
                       [~,EN]=ismember([min(Kites(i,[3,4])),max(Kites(i,[3,4]))],Edges,'rows');
                       plot(PointCoordinates(Kites(i,[3,4]),1),PointCoordinates(Kites(i,[3,4]),2),'Color',[1,0,0],'LineWidth',3);
                       V=[0,-1;1,0]*[diff(PointCoordinates(Kites(i,[3,4]),1));diff(PointCoordinates(Kites(i,[3,4]),2))];
                       V=V/norm(V);
                       text(sum(PointCoordinates(Kites(i,[3,4]),1))/2+V(1)*0.1,sum(PointCoordinates(Kites(i,[3,4]),2))/2+V(2)*0.1,int2str( EdgeOrder(EN)))
                   else
                        plot(PointCoordinates(Kites(i,[3,4]),1),PointCoordinates(Kites(i,[3,4]),2),'Color','k','LineWidth',3);
                   end
                   
                   
                   %Checks if red or black
                   [rows1,~]=find(Kites==Kites(i,4));
                   [rows2,~]=find(Kites==Kites(i,1));
                   AK=rows1(ismember(rows1,rows2));
                   if sum(KiteSeen(AK))~=length(AK)
                       %if red, the priority number is calculated and
                       %plotted
                       [~,EN]=ismember([min(Kites(i,[4,1])),max(Kites(i,[4,1]))],Edges,'rows');
                       plot(PointCoordinates(Kites(i,[4,1]),1),PointCoordinates(Kites(i,[4,1]),2),'Color',[1,0,0],'LineWidth',3);
                       V=[0,-1;1,0]*[diff(PointCoordinates(Kites(i,[4,1]),1));diff(PointCoordinates(Kites(i,[4,1]),2))];
                       V=V/norm(V);
                       text(sum(PointCoordinates(Kites(i,[4,1]),1))/2+V(1)*0.1,sum(PointCoordinates(Kites(i,[4,1]),2))/2+V(2)*0.1,int2str( EdgeOrder(EN)))
                   else
                       plot(PointCoordinates(Kites(i,[4,1]),1),PointCoordinates(Kites(i,[4,1]),2),'Color','k','LineWidth',3);
                   end

                   %plots the four edgepoints (one primal, one dual and two
                   %edgepoints
                   plot(PointCoordinates(Kites(i,1),1),PointCoordinates(Kites(i,1),2),'.','Color',[0.5,0,0],'Markersize',40);
                   plot(PointCoordinates(Kites(i,[2,4]),1),PointCoordinates(Kites(i,[2,4]),2),'.','Color',[0,0.5,0],'Markersize',40);
                   plot(PointCoordinates(Kites(i,3),1),PointCoordinates(Kites(i,3),2),'.','Color',[0,0,0.5],'Markersize',40);
                   
                end
            end
            pause(0.1)

            %set the title
            title("Quad Graph")

            %if the figure should be printed...
            if PrintImages
                exportgraphics(gca,Prefix+"QuadGraphBuild"+sum(KiteSeen)+Suffix);
            end
        end
    end
end


%plots several graphics
if visualization
    
    figure

    %plots all lines of the kites
    for i = 1:AmountOfKites
        if KiteSeen(i)==1 
           plot(PointCoordinates(Kites(i,[1,2]),1),PointCoordinates(Kites(i,[1,2]),2),'k','LineWidth',2);
           hold on
           plot(PointCoordinates(Kites(i,[2,3]),1),PointCoordinates(Kites(i,[2,3]),2),'k','LineWidth',2);
           plot(PointCoordinates(Kites(i,[3,4]),1),PointCoordinates(Kites(i,[3,4]),2),'k','LineWidth',2);
           plot(PointCoordinates(Kites(i,[4,1]),1),PointCoordinates(Kites(i,[4,1]),2),'k','LineWidth',2);
        end
    end
    axis equal
    axis off

    %Plots the primal dual and edge points
    plot(PointCoordinates(1:AmountOfVertices,1),PointCoordinates(1:AmountOfVertices,2),'.','Color',[0.5,0,0],'Markersize',20)
    plot(PointCoordinates(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,1),PointCoordinates(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,2),'.','Color',[0,0,0.5],'Markersize',20)
    plot(PointCoordinates(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),'.','Color',[0,0.5,0],'Markersize',20)
    
    %set the title
    title("Quad Graph")

    %if the graphics should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"QuadGraphWithPointsNew"+Suffix);
    end


    %Plots the circles of the circle packing in 2D
    figure

    %plots the primal circles
    for i=1:length(Kites)
        PrimalPoint=PointCoordinates(Kites(i,1,:),:);
        CirclePoint=PrimalPoint;
        EdgePoint=PointCoordinates(Kites(i,2,:),:);
        Radius=norm(CirclePoint-EdgePoint);
        Circle=computeCircle3D([0,0,1],Radius,[CirclePoint,0],50);
        plot(Circle(:,1),Circle(:,2),':','LineWidth',2,'Color',[0.5,0,0])
        hold on
    end

    %plots the dual circles
    for i=1:length(Kites)
        PrimalPoint=PointCoordinates(Kites(i,3,:),:);
        CirclePoint=PrimalPoint;
        EdgePoint=PointCoordinates(Kites(i,4,:),:);
        Radius=norm(CirclePoint-EdgePoint);
        Circle=computeCircle3D([0,0,1],Radius,[CirclePoint,0],50);
        plot(Circle(:,1),Circle(:,2),'-','LineWidth',2,'Color',[0.0,0,0.5])
        hold on
    end
    



    %plots the four straight line circles
    minx=min(PointCoordinates(:,1));
    miny=min(PointCoordinates(:,2));

    maxx=max(PointCoordinates(:,1));
    maxy=max(PointCoordinates(:,2));


    plot([maxx,maxx],[miny-1,maxy+1],'-','LineWidth',2,'Color',[0.0,0,0.5])
    plot([minx,minx],[miny-1,maxy+1],'-','LineWidth',2,'Color',[0.0,0,0.5])


    plot([minx-1,maxx+1],[maxy,maxy],':','LineWidth',2,'Color',[0.5,0,0])
    plot([minx-1,maxx+1],[miny,miny],':','LineWidth',2,'Color',[0.5,0,0])


    axis equal
    axis off

    %set the title
    title("2D Circle Packing")

    %if the plot should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"Circle2D"+Suffix);
    end


end

%Status that Kites are drawn
if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Kites were drawn.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end


%-------------------------------------------------------------------------
%       Compute the inverse stereographic projection -> Algorithm 8 of
%       Dissertation
%-------------------------------------------------------------------------

%We start by computing the point distance for each point.
%The point distance is the average distance form the point to its
%neighbored points

%Initiation the point distance by the radii (For the primal and dual points
%it is that radius. For the edge point we have to compute the distance).
PointDistance=RPoints;

%The amount of neighbored points 
PointDivider=0*RPoints;

%Goes through all kites and add the distance of the edge points. This
%yields to a double counting for each inner line.
for i=1:AmountOfKites
    PointDistance(Kites(i,2))=PointDistance(Kites(i,2))+RPoints(Kites(i,1));
    PointDistance(Kites(i,2))=PointDistance(Kites(i,2))+RPoints(Kites(i,3));
    
    PointDistance(Kites(i,4))=PointDistance(Kites(i,4))+RPoints(Kites(i,1));
    PointDistance(Kites(i,4))=PointDistance(Kites(i,4))+RPoints(Kites(i,3));

    PointDivider(Kites(i,2))=PointDivider(Kites(i,2))+2;
    PointDivider(Kites(i,4))=PointDivider(Kites(i,4))+2;
end

%The non mentioned points are the points which are not drawn. They get the
%maximum value
PointDistance(PointDistance==0)=max(PointDistance);

%Compute the final point Distance for the edge points. 
PointDividerLocal=PointDivider(AmountOfVertices+1:AmountOfVertices+AmountOfEdges);
PointDistanceLocal=PointDistance(AmountOfVertices+1:AmountOfVertices+AmountOfEdges);

%The distance of the inner regular edge points are computed
PointDistanceLocal(PointDividerLocal==8)=PointDistanceLocal(PointDividerLocal==8)/8;

%The distance for points on the boundary or non drawn point is set to the
%maximum
PointDistanceLocal(PointDividerLocal<8)=max(PointDistance);

%The values are set to the big list
PointDistance(AmountOfVertices+1:AmountOfVertices+AmountOfEdges)=PointDistanceLocal;

%-------------------------------------------------------------------------
%Creates a list of possible starting points. These list just consists out
%of drawn points

%The coordinates
PointCoordinatesStartValue=PointCoordinates(logical(PointSeen),:);

%The point distance
PointDistanceStartValue=PointDistance(logical(PointSeen));

%The vector for the possible alpha value
PointAlphaStartValue=0*PointDistanceStartValue;

%Creates the list of the points, which will be mapped to the unit sphere
EdgePointCoordinates=PointCoordinates(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,:);
EdgePointSeen=PointSeen(AmountOfVertices+1:AmountOfVertices+AmountOfEdges);

%The final list
EdgePointCoordinatesSmall=EdgePointCoordinates(logical(EdgePointSeen),:);

if status
    clc 
    fprintf(StatusString);
    disp(['(', num2str(statusCounter) , ') Compute the inverse stereographic projection:'])
    disp('    Compute the list of starting points.')
end


%Orders the start points by its point distance
[~,PointPriorityList]=sort(PointDistanceStartValue);


if visualization

    %Plots the range of possible starting points for the algorithm of the
    %inverse stereographic projection
    figure

    %plots the quad graph
    for i = 1:AmountOfKites
        if KiteSeen(i)==1 
           plot(PointCoordinates(Kites(i,[1,2]),1),PointCoordinates(Kites(i,[1,2]),2),'k','LineWidth',2);
           hold on
           plot(PointCoordinates(Kites(i,[2,3]),1),PointCoordinates(Kites(i,[2,3]),2),'k','LineWidth',2);
           plot(PointCoordinates(Kites(i,[3,4]),1),PointCoordinates(Kites(i,[3,4]),2),'k','LineWidth',2);
           plot(PointCoordinates(Kites(i,[4,1]),1),PointCoordinates(Kites(i,[4,1]),2),'k','LineWidth',2);
        end
    end
    axis equal
    axis off
    drawnow

    %goes through a 100x100 grid of the quad graph and tests if for this
    %initial data it is possible to get the parameter for the inverse
    %stereographich projection
    for i = min(PointCoordinates(:,1)):(max(PointCoordinates(:,1))-min(PointCoordinates(:,1)))/100:max(PointCoordinates(:,1))
        for j = min(PointCoordinates(:,2)):(max(PointCoordinates(:,2))-min(PointCoordinates(:,2)))/100:max(PointCoordinates(:,2))
            
            %Sets the initial Vx and Vy values
            CoordinateInfluence=[];
            CoordinateInfluence([1,2])=-[i,j];
        
            %find a right border for the alpha value (0 every point is
            %projected to the southpole, infinity, almost all points are
            %projected to the northpole
            positiveValueFound=false;
            potenz=0;
            leftBorder=0;
            while ~positiveValueFound
                sign=ComputeZCoordinateBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),10^potenz);
                
                %if the z value of the barycenter is greater then 0, a right boarder is found 
                if sign>0
                    rightBorder=10^potenz;
                    positiveValueFound=true;
                else
                    potenz=potenz+1;
                end
            end
        
            %uses the interval bisection method to fix the alpha value to a good
            %starting point
            IntervalMean=ComputeZCoordinateBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),(rightBorder+leftBorder)/2);      
            while max(abs(IntervalMean)) > Tolerance 
                if IntervalMean<0
                    leftBorder=(rightBorder+leftBorder)/2;
                else
                    rightBorder=(rightBorder+leftBorder)/2;
                end
                
                IntervalMean=ComputeZCoordinateBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),(rightBorder+leftBorder)/2);
            end

            %sets the alpha value
            CoordinateInfluence(3)=(rightBorder+leftBorder)/2;
        
            CoordinateInfluence=CoordinateInfluence';
            CoordinateInfluenceStart=CoordinateInfluence;

            %Compute the error of the restrictions
            Gradient=ComputeBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));

            %Set the iterations of Newtons Method (to avoid infinity loop)
            Iterations=1;
            
            %As long as at least one restriction is hurt
            while max(abs(Gradient)) > Tolerance && Iterations < MaxIterations
                
                GradientLast=Gradient;
                            
                %Evaluate Hesse Matrix
                Hesse=ComputeHesseInverseStereographicProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));
            
                %Compute direction (as the system might be linear dependend because
                %of redundand restrictions, warnings are off
                warning off
                Direction=Hesse\(Gradient);
                warning on
                
                %Set new evaluation point
                CoordinateInfluence=CoordinateInfluence-Direction;
            
                %Evaluate new gradient
                Gradient=ComputeBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));
                Iterations=Iterations+1;
                
            end
    
            %checks, if the parameters could be computed with newtons
            %method depending on the given starting point

            %checks, if The barycenter has a nan value
            if sum(isnan(Gradient))>0

                %if so and the last coordinate of the barycenter is greater
                %then 0
                if GradientLast(3)>0
                    %the color is red
                    plot(-CoordinateInfluenceStart(1),-CoordinateInfluenceStart(2),'.r','Markersize',4);
                else
                    %if not the color is blue
                    plot(-CoordinateInfluenceStart(1),-CoordinateInfluenceStart(2),'.b','Markersize',4);

                end

            %Checks, if a value of the barycenter is infinity
            elseif sum(isinf(Gradient)) > 0
                %if so, the color is blue
                plot(-CoordinateInfluenceStart(1),-CoordinateInfluenceStart(2),'.y','Markersize',4);

            %Checks if the demanded tolerance is reached
            elseif max(abs(Gradient)) >= Tolerance

                %if not, the color is yellor
                plot(-CoordinateInfluenceStart(1),-CoordinateInfluenceStart(2),'.y','Markersize',4);

            else

                %if so, the mathod was sucessfull and the color is green
                plot(-CoordinateInfluenceStart(1),-CoordinateInfluenceStart(2),'.g','Markersize',4);

            end
        
        end
        drawnow;
    end

    %set the title
    title("Possible Starting Points")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"ConvergenceRadius"+Suffix);
    end

end

%Start values for the while loop
SolutionFound=false;
counter=1;

%goes trough all drawn points and tries to find a solution for the inverse
%stereographic projection
while ~SolutionFound && counter< length(PointPriorityList)
    if status
        clc 
        fprintf(StatusString);
        disp(['(', num2str(statusCounter) , ') Compute the inverse stereographic projection:'])
        disp('    Try to find the values for Vx, Vy, alpha by Newtons method.' )
        disp(['    Point: ' , num2str(counter) ,' / ' , num2str(length(PointPriorityList)) ,'.'])
    end

    %Set the starting values for Vx and Vy
    CoordinateInfluence=[];
    CoordinateInfluence([1,2])=-PointCoordinatesStartValue(PointPriorityList(counter,:),:);

    %One start value for the interval bisection method
    leftBorder=0;

    %Computes  the second value for for the inteval
    %bisection method a value for alpha that is greater than zero 
    positiveValueFound=false;
    rightBorder=1;
    
    %-------------------------------------------------------------------
    %Algorithm 9 of Dissertation

    %multiplyes the right border by 10 as long as the z value of the barycenter
    %is greater than zero.
    while ~positiveValueFound
        sign=ComputeZCoordinateBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),rightBorder);
        if sign>0
            positiveValueFound=true;
        else
            rightBorder=rightBorder*10;
        end
    end

    %Starting now th interval bisection. Computing the z value of the
    %barycenter for the mean value of the 2 borders of alpha 
    IntervalMean=ComputeZCoordinateBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),(rightBorder+leftBorder)/2);
    
    %Do the interval bisection as long as the third coordinate of the
    %barycenter is smaller than the given tolerance
    while abs(IntervalMean) > Tolerance 

        %rearanging the borders
        if IntervalMean<0
            leftBorder=(rightBorder+leftBorder)/2;
        else
            rightBorder=(rightBorder+leftBorder)/2;
        end

        %compute the new z value
        IntervalMean=ComputeZCoordinateBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),(rightBorder+leftBorder)/2);
    end

    %-------------------------------------------------------------------

    %Set the start value of alpha to the list (in case we need it for the
    %fsolve command
    PointAlphaStartValue(PointPriorityList(counter))=(rightBorder+leftBorder)/2;

    %Set the alpha value for Newton's method
    CoordinateInfluence(3)=PointAlphaStartValue(PointPriorityList(counter));
    CoordinateInfluence=CoordinateInfluence';

    %Compute the barycenter of the inverse sterreographic projection
    Gradient=ComputeBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));
     
    %Set the iterations of Newtons Method (to avoid infinity loops)
    Iterations=1;
    
    %As long as at least one restriction is hurt
    while max(abs(Gradient)) >= Tolerance && Iterations < MaxIterations && sum(isnan(Gradient))==0 && sum(isnan(Gradient))==0
        
        if status
            clc 
            fprintf(StatusString);
            disp(['(', num2str(statusCounter) , ') Compute the inverse stereographic projection:'])
            disp('    Try to find the values for Vx, Vy, alpha by Newtons method.' )
            disp(['    Point: ' , num2str(counter) ,' / ' , num2str(length(PointPriorityList)) ,'.'])
            disp(['    Iteration: ', num2str(Iterations),' / ',num2str(MaxIterations)  ])
            disp(['    Error: ', num2str(max(abs(Gradient))),' / ',num2str(Tolerance) ])
        end
            
        %Evaluate Hesse matrix
        Hesse=ComputeHesseInverseStereographicProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));
    
    
        %Compute direction (as the system might be linear dependend because
        %of redundand restrictions, warnings are off
        warning off
        Direction=Hesse\(Gradient);
        warning on
        
        %Set new evaluation point
        CoordinateInfluence=CoordinateInfluence-Direction;
    
        %Evaluate new gradient
        Gradient=ComputeBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));
        Iterations=Iterations+1;
        
    end
    
    %If a solution was found we leave the while loop. If not we try another start point 
    if max(abs(Gradient))<Tolerance && sum(isnan(Gradient))==0 && sum(isnan(Gradient))==0
        SolutionFound=true;
    end

    counter=counter+1;
end

%If there was no solution found with Newton's method, try to find a
%solution with fsolve
if SolutionFound == false
    
    counter=1;
    
    %Set the options for the standard solver
    options = optimoptions('fsolve','FunctionTolerance',Tolerance,'StepTolerance',Tolerance,'OptimalityTolerance',Tolerance,'Display','off');
    
    %goes trough all drawn points and tries to find a solution for the inverse
    %stereographic projection
    while ~SolutionFound && counter< length(PointPriorityList)
        if status
            clc 
            fprintf(StatusString);
            disp(['(', num2str(statusCounter) , ') Compute the inverse stereographic projection:'])
            disp('    Newtons method failed for all starting points.' )
            disp('    Try to find the values for Vx, Vy, alpha by Matlabs fsolve.' )
            disp(['    Point: ' , num2str(counter) ,' / ' , num2str(length(PointPriorityList)) ,'.'])
        end
        
        %Set the starting values for Vx and Vy
        CoordinateInfluence=[];
        CoordinateInfluence([1,2])=-PointCoordinatesStartValue(PointPriorityList(counter,:),:);
    
        %Set the alpha value
        CoordinateInfluence(3)=PointAlphaStartValue(PointPriorityList(counter));

    
        %Define function for fsolve
        f=@(x) ComputeBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),x(1),x(2),x(3));

        %Use matlab standard solver to solve the system 
        [CoordinateInfluence,~,~]=fsolve(f,CoordinateInfluence,options);

        %compute the barycenter
        Gradient=ComputeBarycenterInverseProjection(EdgePointCoordinatesSmall(:,1),EdgePointCoordinatesSmall(:,2),...
                CoordinateInfluence(1),CoordinateInfluence(2),CoordinateInfluence(3));
        
        %If a solution was found we leave the while loop. If not we try another start point 
        if max(abs(Gradient))<Tolerance && sum(isnan(Gradient))==0 && sum(isnan(Gradient))==0
            SolutionFound=true;
        end
    
        counter=counter+1;
    end

end

%if still no solution was found
if ~SolutionFound
    error('It was not possible to compute the inverse stereographic projection.')
end

%Status that the inverse stereographic projection was found
if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') Inverse stereographic projection was found.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end



%-------------------------------------------------------------------------
%       Compute the primal and dual polygon
%-------------------------------------------------------------------------

%The 3D coordinates of all points (primal edge and dual)
PointCoordinates3D=zeros(AmountOfVertices+AmountOfEdges+AmountOfFaces,3);

%Set the north pole (the first edgepoint)
PointCoordinates3D(AmountOfVertices+1,:)=[0,0,1];

%Set the values for the inverse stereographic projection
X=PointCoordinates(AmountOfVertices+2:AmountOfVertices+AmountOfEdges,1);
Y=PointCoordinates(AmountOfVertices+2:AmountOfVertices+AmountOfEdges,2);
Vx=CoordinateInfluence(1);
Vy=CoordinateInfluence(2);
alpha=CoordinateInfluence(3);



%Compute the values for the edge points
x=(2*((X+Vx)*alpha))./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);
y=(2*((Y+Vy)*alpha))./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);
z=(-1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2)./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);

%Set the edge points
PointCoordinates3D(AmountOfVertices+2:AmountOfVertices+AmountOfEdges,:)=[x,y,z];

%Compute the primal points by solving a linear system consisting out of
%planes
for i = 1:AmountOfVertices
    [ActiveKites,~]=find(KiteAll==i);
    ActiveEdgePoints=unique(KiteAll(ActiveKites,[2,4]));
    PointGLS=PointCoordinates3D(ActiveEdgePoints,:);
    PointRightSide=ones(length(ActiveKites),1);

    %If the system is 3x3 we dont need least squares
    if length(ActiveKites)==3
        PointCoordinates3D(i,:)=(PointGLS\PointRightSide)';
    else
        PointCoordinates3D(i,:)=((PointGLS'*PointGLS)\(PointGLS'*PointRightSide))';
    end
end

%Compute the dual points by solving a linear system consisting out of
%planes
for i = AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces
    [ActiveKites,~]=find(KiteAll==i);
    ActiveEdgePoints=unique(KiteAll(ActiveKites,[2,4]));
    PointGLS=PointCoordinates3D(ActiveEdgePoints,:);
    PointRightSide=ones(length(ActiveKites),1);

        %If the system is 3x3 we dont need least squares
    if length(ActiveKites)==3
        PointCoordinates3D(i,:)=(PointGLS\PointRightSide)';
    else
        PointCoordinates3D(i,:)=((PointGLS'*PointGLS)\(PointGLS'*PointRightSide))';
    end
end

%plots the quad graph with the center and the radius of the inverse
%stereographic projection
if visualization

    figure

    %moves and scales the points
    PointCoordinatesLocal=(PointCoordinates+[Vx,Vy])*alpha;

    %plot all the kites
    for i = 1:AmountOfKites
        if KiteSeen(i)==1 
           plot(PointCoordinatesLocal(Kites(i,[1,2]),1),PointCoordinatesLocal(Kites(i,[1,2]),2),'k','LineWidth',2);
           hold on
           plot(PointCoordinatesLocal(Kites(i,[2,3]),1),PointCoordinatesLocal(Kites(i,[2,3]),2),'k','LineWidth',2);
           plot(PointCoordinatesLocal(Kites(i,[3,4]),1),PointCoordinatesLocal(Kites(i,[3,4]),2),'k','LineWidth',2);
           plot(PointCoordinatesLocal(Kites(i,[4,1]),1),PointCoordinatesLocal(Kites(i,[4,1]),2),'k','LineWidth',2);
        end
    end
    axis equal

    %plot the center after movement
    plot(0,0,'.','Markersize',20,'Color',[0.5,0.5,1])

    %plot the radius (all points on the radius are mapped to the equator)
    rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',3,'EdgeColor',[0.5,0.5,1]);
    axis off

    %set the title
    title("Quad Graph with Center and Radius")

    %if graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"QuadGraphAdjusted"+Suffix);
    end

    %also plot the edge points
    plot(PointCoordinatesLocal(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinatesLocal(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),'.','Markersize',20,'Color',[0,0.5,0])
    
    %replot center and radius (prevent overlapping of edgepoints)
    plot(0,0,'.','Markersize',40,'Color',[0.5,0.5,1])
    rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',3,'EdgeColor',[0.5,0.5,1]);
    
    %set the title
    title("Quad Graph with Center and Radius")

    %if graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"QuadGraphWithPointsAdjusted"+Suffix);
    end
end

%plots the quad graph and the projected points to the sphere
if visualization
    figure

    %shift local points
    PointCoordinatesLocal=(PointCoordinates+[Vx,Vy])*alpha;
    for i = 1:AmountOfKites
        if KiteSeen(i)==1 

           %plot the kites
           plot(PointCoordinatesLocal(Kites(i,[1,2]),1),PointCoordinatesLocal(Kites(i,[1,2]),2),'Color',[0,0,0],'LineWidth',2);
           hold on
           plot(PointCoordinatesLocal(Kites(i,[2,3]),1),PointCoordinatesLocal(Kites(i,[2,3]),2),'Color',[0,0,0],'LineWidth',2);
           plot(PointCoordinatesLocal(Kites(i,[3,4]),1),PointCoordinatesLocal(Kites(i,[3,4]),2),'Color',[0,0,0],'LineWidth',2);
           plot(PointCoordinatesLocal(Kites(i,[4,1]),1),PointCoordinatesLocal(Kites(i,[4,1]),2),'Color',[0,0,0],'LineWidth',2);


           %plot the edges from  the quad graph to the northpole
           plot3([PointCoordinatesLocal(Kites(i,2),1),0],[PointCoordinatesLocal(Kites(i,2),2),0],[0,1],'Color',[0,0.5,0],'LineWidth',1);
           plot3([PointCoordinatesLocal(Kites(i,4),1),0],[PointCoordinatesLocal(Kites(i,4),2),0],[0,1],'Color',[0,0.5,0],'LineWidth',1);

           %plot the edges from the projected points to the northpole
           plot3([PointCoordinates3D(Kites(i,2),1),0],[PointCoordinates3D(Kites(i,2),2),0],[PointCoordinates3D(Kites(i,2),3),1],'Color',[0,0.5,0],'LineWidth',1);
           plot3([PointCoordinates3D(Kites(i,4),1),0],[PointCoordinates3D(Kites(i,4),2),0],[PointCoordinates3D(Kites(i,4),3),1],'Color',[0,0.5,0],'LineWidth',1);
            
        end
    end
    axis equal
    axis off

    %plot the edge points of the quad graph
    plot(PointCoordinatesLocal(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinatesLocal(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),'.','Color',[0,0.5,0],'Markersize',25)
    
    %plot the projected points on the sphere
    plot3(PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize',10);
    
    %adjust the point of view
    view(View);
    [X,Y,Z]=sphere;
    CO(:,:,1) = zeros(21)+0.5; % red
    CO(:,:,2) = zeros(21)+0.5; % green
    CO(:,:,3) = ones(21).*1; % blue

    %plots the sphere
    surf(X,Y,Z,CO,'FaceAlpha','0.2')

    %set the title
    title("Projection to Sphere")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"ProjectionToSphere"+Suffix);
    end

end

%plots the edge points on the sphere
if visualization
    figure

    %plot the edge points
    plot3(PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize',25);
    hold on

    %adjust the view
    view(View);

    %plot the sphere
    [X,Y,Z]=sphere;
    CO(:,:,1) = zeros(21)+0.5; % red
    CO(:,:,2) = zeros(21)+0.5; % green
    CO(:,:,3) = ones(21).*1; % blue
    surf(X,Y,Z,CO,'FaceAlpha','0.2')
   
    axis equal
    axis off

    %set the title
    title("Projected Points")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"EdgePointSphere"+Suffix);
    end

end


% plot the primal polytope including the edgepoints
if visualization

    figure
    
    %plot the lines of the primal polytope
    for i = 1:length(KiteAll)
           plot3(PointCoordinates3D(KiteAll(i,[1,2]),1),PointCoordinates3D(KiteAll(i,[1,2]),2),PointCoordinates3D(KiteAll(i,[1,2]),3),'Color',[0.5,0,0],'LineWidth',2);
           hold on
      
           plot3(PointCoordinates3D(KiteAll(i,[4,1]),1),PointCoordinates3D(KiteAll(i,[4,1]),2),PointCoordinates3D(KiteAll(i,[4,1]),3),'Color',[0.5,0,0],'LineWidth',2);
    end
    axis equal
    axis off

    %plot the edgepoints
    plot3(PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize',25);
    
    %plot the primal points
    plot3(PointCoordinates3D(1:AmountOfVertices,1),PointCoordinates3D(1:AmountOfVertices,2),PointCoordinates3D(1:AmountOfVertices,3),'.','Color',[0.5,0,0],'Markersize',25);
    
    %adjust the view
    view(View);

    %set the title
    title("Primal Polytope")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"PolytopePrimal"+Suffix);
    end
    

end


%plot the dual polytope and the edgepoints
if visualization
    
    figure
    
    %plot the lines of the polytope
    for i = 1:length(KiteAll)
       plot3(PointCoordinates3D(KiteAll(i,[2,3]),1),PointCoordinates3D(KiteAll(i,[2,3]),2),PointCoordinates3D(KiteAll(i,[2,3]),3),'Color',[0,0,0.5],'LineWidth',2);
       hold on
       plot3(PointCoordinates3D(KiteAll(i,[3,4]),1),PointCoordinates3D(KiteAll(i,[3,4]),2),PointCoordinates3D(KiteAll(i,[3,4]),3),'Color',[0,0,0.5],'LineWidth',2);
    end
    axis equal
    axis off

    %plot the edge points
    plot3(PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize', 25);
    
    %plot the dual points
    plot3(PointCoordinates3D(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,1),PointCoordinates3D(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,2),PointCoordinates3D(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,3),'.','Color',[0,0,0.5],'Markersize',25);

    %adjust the view
    view(View);

    %set the title
    title("Dual Polytope")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"PolytopeDual"+Suffix);
    end
end



%plot the primal and dual polytope in one image
if visualization
    figure
    
    %plot the edges of the primal and dual polytope
    for i = 1:length(KiteAll)
           plot3(PointCoordinates3D(KiteAll(i,[1,2]),1),PointCoordinates3D(KiteAll(i,[1,2]),2),PointCoordinates3D(KiteAll(i,[1,2]),3),'LineWidth',2,'Color',[0.5,0,0]);
           hold on
           plot3(PointCoordinates3D(KiteAll(i,[2,3]),1),PointCoordinates3D(KiteAll(i,[2,3]),2),PointCoordinates3D(KiteAll(i,[2,3]),3),'LineWidth',2,'Color',[0,0,0.5]);
           plot3(PointCoordinates3D(KiteAll(i,[3,4]),1),PointCoordinates3D(KiteAll(i,[3,4]),2),PointCoordinates3D(KiteAll(i,[3,4]),3),'LineWidth',2,'Color',[0,0,0.5]);
           plot3(PointCoordinates3D(KiteAll(i,[4,1]),1),PointCoordinates3D(KiteAll(i,[4,1]),2),PointCoordinates3D(KiteAll(i,[4,1]),3),'LineWidth',2,'Color',[0.5,0,0]);
    end
    axis equal
    axis off

    %plot the edge points
    plot3(PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3D(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize',25);
    
    %plot the primal points
    plot3(PointCoordinates3D(1:AmountOfVertices,1),PointCoordinates3D(1:AmountOfVertices,2),PointCoordinates3D(1:AmountOfVertices,3),'.','Color',[0.5,0,0],'Markersize',25);
    
    %plot the dual points
    plot3(PointCoordinates3D(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,1),PointCoordinates3D(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,2),PointCoordinates3D(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,3),'.','Color',[0,0,0.5],'Markersize',25);
    
    %adjust the view
    view(View);

    %set the title
    title("Primal and Dual Polytope")

    %if the graphic schould be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"PolytopeBoth"+Suffix);
    end

end


%plot the primal polytope and the dual points projected on the faces of the
%primal polytope
if visualization
    figure

    PointCoordinates3DTemp=PointCoordinates3D;

    %projects the dual points
    PointCoordinates3DTemp(KiteAll(:,3),:)=1./(PointCoordinates3DTemp(KiteAll(:,3),1).^2+PointCoordinates3DTemp(KiteAll(:,3),2).^2+PointCoordinates3DTemp(KiteAll(:,3),3).^2).*PointCoordinates3DTemp(KiteAll(:,3),:);
    
    %plot the lines of the primal and dual polytope
    for i = 1:length(KiteAll)
           plot3(PointCoordinates3DTemp(KiteAll(i,[1,2]),1),PointCoordinates3DTemp(KiteAll(i,[1,2]),2),PointCoordinates3DTemp(KiteAll(i,[1,2]),3),'LineWidth',2,'Color',[0.5,0,0]);
           hold on
           plot3(PointCoordinates3DTemp(KiteAll(i,[2,3]),1),PointCoordinates3DTemp(KiteAll(i,[2,3]),2),PointCoordinates3DTemp(KiteAll(i,[2,3]),3),'LineWidth',2,'Color',[0,0,0.5]);
           plot3(PointCoordinates3DTemp(KiteAll(i,[3,4]),1),PointCoordinates3DTemp(KiteAll(i,[3,4]),2),PointCoordinates3DTemp(KiteAll(i,[3,4]),3),'LineWidth',2,'Color',[0,0,0.5]);
           plot3(PointCoordinates3DTemp(KiteAll(i,[4,1]),1),PointCoordinates3DTemp(KiteAll(i,[4,1]),2),PointCoordinates3DTemp(KiteAll(i,[4,1]),3),'LineWidth',2,'Color',[0.5,0,0]);
    end
    axis equal
    axis off

    %plot the edgepoints
    plot3(PointCoordinates3DTemp(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,1),PointCoordinates3DTemp(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,2),PointCoordinates3DTemp(AmountOfVertices+1:AmountOfVertices+AmountOfEdges,3),'.','Color',[0,0.5,0],'Markersize',25);
    
    %plot the primal points
    plot3(PointCoordinates3DTemp(1:AmountOfVertices,1),PointCoordinates3DTemp(1:AmountOfVertices,2),PointCoordinates3DTemp(1:AmountOfVertices,3),'.','Color',[0.5,0,0],'Markersize',25);
    
    %plot the projected dual points
    plot3(PointCoordinates3DTemp(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,1),PointCoordinates3DTemp(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,2),PointCoordinates3DTemp(AmountOfVertices+AmountOfEdges+1:AmountOfVertices+AmountOfEdges+AmountOfFaces,3),'.','Color',[0,0,0.5],'Markersize',25);
    
    %adjust the view
    view(View);

    %set the title
    title("Primal Polytope with Projected Dual Polytope")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"PolytopeBothProjected"+Suffix);
    end
    
end


%plot the circle packing in 3D
if visualization

    PointCoordinates3DTemp=PointCoordinates3D;
    figure

    %plot the primal circles
    for i=1:AmountOfVertices

        %get the current primal point
        PrimalPoint=PointCoordinates3DTemp(i,:);

        %get the projection if this point, which is the center of the
        %circle
        CirclePoint=PrimalPoint/norm(PrimalPoint)^2;

        %gets an edgepoint to compute the radius
        [row,col]=find(KiteAll==i);
        row=row(1);
        col=col(1);
        if col==4
            EdgePointIndex=KiteAll(row,3);
        else
            EdgePointIndex=KiteAll(row,col+1);
        end
        EdgePoint=PointCoordinates3DTemp(EdgePointIndex,:);

        %compute the radius
        Radius=norm(CirclePoint-EdgePoint);

        %compute the circle
        Circle=computeCircle3D(CirclePoint,Radius,CirclePoint,50);

        %plot the circle
        plot3(Circle(:,1),Circle(:,2),Circle(:,3),':','LineWidth',2,'Color',[0.5,0,0])
        hold on
    end

    %plot the dual circles
    for i=AmountOfVertices+AmountOfEdges+1:length(PointCoordinates3DTemp)

        %get the current dual point
        DualPoint=PointCoordinates3DTemp(i,:);

        %get its projection, which is the center of the circle
        CirclePoint=DualPoint/norm(DualPoint)^2;

        %find an edgepoint to compute the radius
        [row,col]=find(KiteAll==i);
        row=row(1);
        col=col(1);
        if col==4
            EdgePointIndex=KiteAll(row,3);
        else
            EdgePointIndex=KiteAll(row,col+1);
        end
        EdgePoint=PointCoordinates3DTemp(EdgePointIndex,:);

        %compute the radius
        Radius=norm(CirclePoint-EdgePoint);

        %compute the circle
        Circle=computeCircle3D(CirclePoint,Radius,CirclePoint,50);

        %plot the circle
        plot3(Circle(:,1),Circle(:,2),Circle(:,3),'-','LineWidth',2,'Color',[0.0,0,0.5])
    end
    
    axis equal
    axis off
    
    %adjust the view
    view(View);

    %set the title
    title("3D Circle Packing")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"Circle3D"+Suffix);
    end

end

%plot only the primal polytope
if visualization
    
    figure
    
    %plot the lines of the primal polytope
    for i = 1:length(KiteAll)
           plot3(PointCoordinates3D(KiteAll(i,[1,2]),1),PointCoordinates3D(KiteAll(i,[1,2]),2),PointCoordinates3D(KiteAll(i,[1,2]),3),'LineWidth',2,'Color',[0.2,0.2,0.2]);
           hold on
           plot3(PointCoordinates3D(KiteAll(i,[4,1]),1),PointCoordinates3D(KiteAll(i,[4,1]),2),PointCoordinates3D(KiteAll(i,[4,1]),3),'LineWidth',2,'Color',[0.2,0.2,0.2]);
    end
    axis equal
    axis off
   
    %plot the primal points
    plot3(PointCoordinates3D(1:AmountOfVertices,1),PointCoordinates3D(1:AmountOfVertices,2),PointCoordinates3D(1:AmountOfVertices,3),'.r','Markersize',40,'Color',[0.5,0.0,0.0]);

    %adjust the view
    view(View);

    %set the title
    title("Primal Polytope")

    %if the graphic should be printed...
    if PrintImages
        exportgraphics(gca,Prefix+"PolytopePrimalColor"+Suffix);
    end


end

if status
    clc 
    StatusString=StatusString+['(', num2str(statusCounter) , ') 3D-Coordinates were computed.\n'];
    fprintf(StatusString);
    statusCounter=statusCounter+1;
end



end


function [V] = ComputeBarycenterInverseProjection(X,Y,Vx,Vy,alpha)
%ComputeBarycenterInverseProjection Evaluates the barycenter of the inverse 
%stereographic projection.
%
%Input: x-coordinates of the points: X
%       y-coordinates of the points: Y
%       Shift in x-direction: Vx
%       Shift in y-direction: Vy
%       Scaling value: alpha       
%
%Output:    The barycenter of the inverse stereographic projection: V

x=[(2*((X+Vx)*alpha))./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);0];
y=[(2*((Y+Vy)*alpha))./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);0];
z=[(-1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2)./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);1];

V=[x,y,z];
V=sum(V)/length(V);
V=V';
end


function [P] = computeCircle3D(NormalVector,Radius,Shift,AmountOfPoints)
%COMPUTECIRCLE3D computes the points of a circle with midpoint Shift,
%radius Radius and normalvector Normalvector.
%
%Input: normalvector: NormalVector
%       radius: Radius
%       midpoint: Shift
%       amount of points in the circle: AmountOfPoints
%
%Output: P: the points on the circle
%
%computeCircle3D(NormalVector,Radius,Shift,AmountOfPoints) computes the 
%points of a circle with midpoint Shift, radius Radius and normalvector 
%Normalvector.

%computes equidistant points on a circle in R^2 with midpoint [0,0] and
%radius Radius
P=[Radius*cos(2*pi*(1:AmountOfPoints)/AmountOfPoints)',Radius*sin(2*pi*(1:AmountOfPoints)/AmountOfPoints)',0*(1:AmountOfPoints)'];

%compute the angle of rotation for this circle
Vektor1=[0,0,1];
Vektor2=NormalVector;
Angle=acos(dot(Vektor1,Vektor2)/norm(Vektor2));

%computs the normal of ratation for this circle
RN=cross(Vektor2,Vektor1);
if norm(RN)~=0
RN=RN/norm(RN);
end

%computs the rotation matrix
RotationMatrix=[RN(1)^2*(1-cos(Angle))+cos(Angle),                RN(1)*RN(2)*(1-cos(Angle))-RN(3)*sin(Angle),          RN(1)*RN(3)*(1-cos(Angle))+RN(2)*sin(Angle) ;...
                RN(1)*RN(2)*(1-cos(Angle))+RN(3)*sin(Angle)       RN(2)^2*(1-cos(Angle))+cos(Angle)                     RN(2)*RN(3)*(1-cos(Angle))-RN(1)*sin(Angle);...
                RN(1)*RN(3)*(1-cos(Angle))-RN(2)*sin(Angle)       RN(2)*RN(3)*(1-cos(Angle))+RN(1)*sin(Angle)           RN(3)^2*(1-cos(Angle))+cos(Angle)      ];
%rotate the points
P=P*RotationMatrix;

%shift the points
P=P+Shift;

%double the first point to get a closed circle
P=[P;P(1,:)];

end


function [Hesse] = ComputeHesseInverseStereographicProjection(X,Y,Vx,Vy,alpha)
%ComputeHesseInverseStereographicProjection Evaluates the Hesse of 
%the inverse stereographic projection.
%
%Input: x-coordinates of the points: X
%       y-coordinates of the points: Y
%       Shift in x-direction: Vx
%       Shift in y-direction: Vy
%       Scaling value: alpha       
%
%Output:    The Hesse of the inverse stereographic projection: Hesse

%Compute the derivative of the nominator of the x coordinate regarding Vx
%Vy and alpha
z1DVx=0.*X+2*alpha;
z1DVy=0.*X;
z1Dalpha=2*(X+Vx);

%Compute the derivative of the nominator of the y coordinate regarding Vx
%Vy and alpha
z2DVx=0.*X;
z2DVy=0.*X+2*alpha;
z2Dalpha=2*(Y+Vy);

%Compute the derivative of the nominator of the z coordinate regarding Vx
%Vy and alpha
z3DVx=2*alpha^2*X+2*alpha^2*Vx;
z3DVy=2*alpha^2*Y+2*alpha^2*Vy;
z3Dalpha=2*alpha*((X+Vx).^2+(Y+Vy).^2);

%Compute the derivative of the denominator regarding VX
%Vy and alpha
NDVx=2*alpha^2*X+2*alpha^2*Vx;
NDVy=2*alpha^2*Y+2*alpha^2*Vy;
NDalpha=2*alpha*((X+Vx).^2+(Y+Vy).^2);

%Compute the nominators of the x, y and z coordinate
z1=(2*(X+Vx)*alpha);
z2=(2*(Y+Vy)*alpha);
z3=((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2-1;

%Compute the denominator
n=((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2+1;

%Compute the Hesse matrix by the quotient rule
Hesse=zeros(3);

Hesse(1,1)=sum((z1DVx.*n-NDVx.*z1)./(n.^2));
Hesse(1,2)=sum((z1DVy.*n-NDVy.*z1)./(n.^2));
Hesse(1,3)=sum((z1Dalpha.*n-NDalpha.*z1)./(n.^2));

Hesse(2,1)=sum((z2DVx.*n-NDVx.*z2)./(n.^2));
Hesse(2,2)=sum((z2DVy.*n-NDVy.*z2)./(n.^2));
Hesse(2,3)=sum((z2Dalpha.*n-NDalpha.*z2)./(n.^2));

Hesse(3,1)=sum((z3DVx.*n-NDVx.*z3)./(n.^2));
Hesse(3,2)=sum((z3DVy.*n-NDVy.*z3)./(n.^2));
Hesse(3,3)=sum((z3Dalpha.*n-NDalpha.*z3)./(n.^2));

%Divide the sum by the amount of points as it is the Hesse of the
%barycenter
Hesse=Hesse/(length(X)+1);
end



function [z] = ComputeZCoordinateBarycenterInverseProjection(X,Y,Vx,Vy,alpha)
%ComputeZCoordinateBarycenterInverseProjection Evaluates the z coordinate 
% of the barycenter of the inverse stereographic projection.
%
%Input: x-coordinates of the points: X
%       y-coordinates of the points: Y
%       Shift in x-direction: Vx
%       Shift in y-direction: Vy
%       Scaling value: alpha       
%
%Output:    The z coordinate of the barycenter of the 
%           inverse stereographic projection: z

z=[(-1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2)./(1+((X+Vx)*alpha).^2+((Y+Vy)*alpha).^2);1];
z=sum(z)/length(z);
end



function [Evaluation] = EvaluateDistance(PointX,PointY,PointCombination,rr)
%EVALUATEDISTANCE evaluates the error of the real distance between 2 points
%and the distance they should have.
%
%Input: x-coordinates of the points: PointX
%       y-coordinates of the points: PointY
%       pair of points which shares a distance restriction: PointCombination
%       right site resp. the distances 2 points should have: rr
%       
%
%Output: Error between the real distance and the distance the points should have: Evaluation 

    Evaluation=sqrt(  (PointX(PointCombination(:,1))-PointX(PointCombination(:,2))).^2+ (PointY(PointCombination(:,1))-PointY(PointCombination(:,2))).^2  )-rr;
end




function [Evaluation] = EvaluateDistanceDerivative(PointX,PointY,PointCombination)
%EvaluateDistanceDerivatives evaluates the Hesse of the error of the real 
% distance between 2 points and the distance they should have.
%
%Input: x-coordinates of the points: PointX
%       y-coordinates of the points: PointY
%       pair of points which shares a distance restriction: PointCombination
%       
%
%Output: The Hesse Matrix of the original linear system 

%Creates the Hesse matrix  
Evaluation=zeros(length(PointCombination),2*length(PointX));

%For each row, there are 4 entries different form zero. The derivative in
%direction a,b,x and y. 

%All derivatives in direction a
Evaluation((PointCombination(:,1)-1)*length(PointCombination)+(1:length(PointCombination))')=...
    (PointX(PointCombination(:,1))-PointX(PointCombination(:,2)))./...
    sqrt(  (PointX(PointCombination(:,1))-PointX(PointCombination(:,2))).^2+ (PointY(PointCombination(:,1))-PointY(PointCombination(:,2))).^2 );

%All derivatives in direction x
Evaluation((PointCombination(:,2)-1)*length(PointCombination)+(1:length(PointCombination))')=...
    (PointX(PointCombination(:,2))-PointX(PointCombination(:,1)))./...
    sqrt(  (PointX(PointCombination(:,1))-PointX(PointCombination(:,2))).^2+ (PointY(PointCombination(:,1))-PointY(PointCombination(:,2))).^2 );

%All derivatives in direction b
Evaluation((PointCombination(:,1)-1+length(PointX))*length(PointCombination)+(1:length(PointCombination))')=...
    (PointY(PointCombination(:,1))-PointY(PointCombination(:,2)))./...
    sqrt(  (PointX(PointCombination(:,1))-PointX(PointCombination(:,2))).^2+ (PointY(PointCombination(:,1))-PointY(PointCombination(:,2))).^2 );

%All derivatives in direction y
Evaluation((PointCombination(:,2)-1+length(PointX))*length(PointCombination)+(1:length(PointCombination))')=...
    (PointY(PointCombination(:,2))-PointY(PointCombination(:,1)))./...
    sqrt(  (PointX(PointCombination(:,1))-PointX(PointCombination(:,2))).^2+ (PointY(PointCombination(:,1))-PointY(PointCombination(:,2))).^2 );

end



function [Evaluation] = EvaluateSRho(EvaluationPoint,RhoCombinations,Pi,Normalization)
%EVALUATESRHO evaluates the system of rhos at the evaluation point
%EvaluationPoint. 
%
%Input: Evaluation Point: EvaluationPoint
%       List of kite combinations: RhoCombinations
%       Right side of the linear system: Pi
%       Information if the sum should be into the vector or not:
%       Normalization
%
%Output: the evaluation of the linear system: Evaluation
%
%EvaluateSRho(EvaluationPoint,RhoCombinations,Pi,Normalization) evaluates the system of 
% rhos at the evaluation point EvaluationPoint.

%Amount of lines of the system
[AmountOfLines,~]=size(EvaluationPoint);



%Return/Result value (with or without sum of all points)
if Normalization
    Evaluation=zeros(AmountOfLines+1,1);
    Evaluation(end)=sum(EvaluationPoint);
else
    Evaluation=zeros(AmountOfLines,1);
end

%Evaluates every line of the system
for i=1:AmountOfLines

    %Finds every kite which belongs to vertex i
    [rows,cols]=find(RhoCombinations==i);

    %Goes through all these kites
    for k=1:length(rows)

        %computes the corresponding other node
        j=RhoCombinations(rows(k),mod(cols(k),2)+1);

        %Adds the sum value
        Evaluation(i)=Evaluation(i)+2*atan(exp(EvaluationPoint(j)-EvaluationPoint(i)));
    end

    %substract the right side
    Evaluation(i)=Evaluation(i)-Pi(i);
end

end



function [Evaluation] = EvaluateSRhoDerivative(EvaluationPoint,RhoCombinations,Normalization)
%EVALUATESRHODERIVATIVE evaluates the Hesse of the system of rhos at the evaluation point
%EvaluationPoint. 
%
%Input: Evaluation Point: EvaluationPoint
%       List of kite combinations: RhoCombinations
%       Information if the sum should be into the vector or not:
%       Normalization
%
%Output: the evaluation of the Hesse of linear system: Evaluation
%
%EvaluateSRho(EvaluationPoint,RhoCombinations,normalization) evaluates the Hesse of the system of 
% rhos at the evaluation point EvaluationPoint.

%Amount of lines of the system
[AmountOfLines,~]=size(EvaluationPoint);


%Return/Result value
if Normalization
    Evaluation=zeros(AmountOfLines+1,AmountOfLines);
    Evaluation(end,:)=1;
else
    Evaluation=zeros(AmountOfLines);
end

%Evaluates every line of the system
for i=1:AmountOfLines

    %Finds every kite which belongs to vertex i
    [rows,cols]=find(RhoCombinations==i);

    %Goes through all these kites
    for k=1:length(rows)

        %computes the corresponding other node
        jlocal=RhoCombinations(rows(k),mod(cols(k),2)+1);

        %Adds the sum value
        Evaluation(i,i)=Evaluation(i,i)-...
                        2 ...
                        *exp(EvaluationPoint(i)+EvaluationPoint(jlocal))...
                        /(exp(2*EvaluationPoint(jlocal))+exp(2*EvaluationPoint(i)));

         Evaluation(i,jlocal)=2 ...
                        *exp(EvaluationPoint(i)+EvaluationPoint(jlocal))...
                        /(exp(2*EvaluationPoint(jlocal))+exp(2*EvaluationPoint(i)));
        
    end
end

end
