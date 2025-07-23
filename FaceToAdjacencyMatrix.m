function [AdjacencyMatrix] = FaceToAdjacencyMatrix(FaceMatrix)
%FACETOADJACENCYMATRIX computes an adjacency matrix out of a face matrix
%by setting the entry (i,j) of the adjacency matrix to 1 if i and j are row 
%neighbored in the FaceMatrix.
%
%Input: FaceMatrix: A matrix of faces (each row represents a face).
%
%Output: AdjacencyMatrix a n x n matrix with entries out of {0,1}. 
%Entry (i,j) is 1 if there is an edge between i and j and 0 if not.  
%
%FaceToAdjacencyMatrix(FaceMatrix) computes an adjacency matrix out of 
%a face matrix.


%The highest entry in the FaceMatrix is the amount of vertices.
AmountOfVertices=max(max(FaceMatrix));

%Creates the adjacency matrix variable
AdjacencyMatrix=zeros(AmountOfVertices);


[AmountOfFaces,~]=size(FaceMatrix);

%Goes through each face 
for i=1:AmountOfFaces
    
    %Look how many points are at one face
    FaceValence=nnz(FaceMatrix(i,:));

    %Goes through all points and identify the neighbors
    for j=1:FaceValence

        %Sets one point
        InitialPont=FaceMatrix(i,j);

        %and the other one as the right neighbor
        if j==FaceValence
            Neighbor=FaceMatrix(i,1);
        else
            Neighbor=FaceMatrix(i,j+1);
        end

        %set the value in the adjacancy matrix
        AdjacencyMatrix(InitialPont,Neighbor)=1;
        AdjacencyMatrix(Neighbor,InitialPont)=1;
    end 

end

end


