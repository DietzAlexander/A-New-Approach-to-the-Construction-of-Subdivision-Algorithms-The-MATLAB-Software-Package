function [S,Lattice] = computeSubdivisionMatrix(Input,Degree,MatrixSize,varargin)
%COMPUTESUVDIVISIONMATRIX Computes a subdivision matrix for the given
%Input. Determines automatically, which Subdivision algorithm has to be
%used.
%
% Input: A number > 2 or a Adjacency, Edge or Face matrix
% Degree: 2 (or 'quadratic', or 'Quadratic') or 3  (or 'cubic', or 'Cubic')
% MatrixSize: 0 (or 'initial', or 'Initial') or 1  (or 'big', or 'Big')
% varargin: Several input arguments (see manual)
%           for variant 1 mu has to be specified lile
%           computeSubdivisionMatrix(~,'mu',1/4)
%           with 0 le mu <1/2.
%           Has especialle specifeid for the 2D case, even if mus is not
%           supported in the 2D case, to get to the whished algorithm
%
% S: Subdivision matrix
% Lattice: Lattice, if the algorithm has this as an output and [] else.

%Checks, if the Degree is 2 or 3
if (isnumeric(Degree) &&  all(Degree == 2)) || strcmp(Degree,'quadratic') || strcmp(Degree,'Quadratic')
    Degree = 2;
elseif (isnumeric(Degree) &&  all(Degree == 3)) || strcmp(Degree,'cubic') || strcmp(Degree,'Cubic')
    Degree = 3;
else
    error('Degree has to be 2 (or quadratic or Quadratic) or 3 (or cubic or Cubic)');
end

%Checks, if The matrix size is initial or big
if (isnumeric(MatrixSize) &&  all(MatrixSize == 0)) || strcmp(MatrixSize,'initial') || strcmp(MatrixSize,'Initial')
    MatrixSize = 0;
elseif (isnumeric(MatrixSize) &&  all(MatrixSize == 1)) || strcmp(MatrixSize,'big') || strcmp(MatrixSize,'Big')
    MatrixSize = 1;
else
    error('MatrixSize has to be 0 (or initial or Initial) or 1 (or big or Big)');
end


%Boolean for the Variant 1 (true) or 2 (false)
V1=false;

%The amount of additional input parameters  
NumberOfAdditionalInput = nargin - 3;

%Checks, if the given keyword could be found in varargin. If yes, values
%were set and the keywords and the values were deleted from varargin. If
%they cannot be found, they are simply ignored.
i=1;
while i <= NumberOfAdditionalInput

    %Specifies whether the eigenstructure shult be plottet or not
    if ischar(varargin{i}) && strcmp(varargin{i},'mu')
        mu=varargin{i+1};
        varargin(i)=[];
        varargin(i)=[];
        NumberOfAdditionalInput=NumberOfAdditionalInput-2;
        V1=true;
    else
        i=i+1;
    end
end

%Check, if mu is numeric
if V1 && ~isnumeric(mu)
    error('mu is not a numeric value.')
end

%Checks if the input is a number (2D case) or not (3D case)
if max(size(Input))==1
    Dimension=2;
else
    Dimension=3;
end

%Set the lattice value. If the algorithm do not have Lattice as an output,
%it remains [];
Lattice=[];

%Check the 12 cases and use the whished algorithm

%2D
if Dimension == 2 && Degree == 2 && MatrixSize== 0 && V1
    S=computeBiQuadraticSubdivisionMatrixV1(Input);
    if mu~=1/4
        warning ('For the 2D case mu~=1/4 is not supported. Set mu=1/4.')
    end

elseif Dimension == 2 && Degree == 2 && MatrixSize== 0 && ~V1
    S=computeBiQuadraticSubdivisionMatrixV2(Input);

elseif Dimension == 2 && Degree == 2 && MatrixSize== 1 && V1
    S=computeBiQuadraticSubdivisionMatrixV1Big(Input);
    if mu~=1/4
        warning ('For the 2D case mu~=1/4 is not supported. Set mu=1/4.')
    end

elseif Dimension == 2 && Degree == 2 && MatrixSize== 1 && ~V1
    S=computeBiQuadraticSubdivisionMatrixV2Big(Input);

elseif Dimension == 2 && Degree == 3 && MatrixSize== 0 
    S = computeBiCubicSubdivisionMatrix(Input);

elseif Dimension == 2 && Degree == 3 && MatrixSize== 1 
    S = computeBiCubicSubdivisionMatrixBig(Input);


%3D
elseif Dimension == 3 && Degree == 2 && MatrixSize== 0 && V1
    S=computeTriQuadraticSubdivisionMatrixV1(Input,mu,varargin{:});
   

elseif Dimension == 3 && Degree == 2 && MatrixSize== 0 && ~V1
    S=computeTriQuadraticSubdivisionMatrixV2(Input,varargin{:});

elseif Dimension == 3 && Degree == 2 && MatrixSize== 1 && V1
    [S,Lattice]=computeTriQuadraticSubdivisionMatrixV1Big(Input,mu,varargin{:});
    

elseif Dimension == 3 && Degree == 2 && MatrixSize== 1 && ~V1
    [S,Lattice]=computeTriQuadraticSubdivisionMatrixV2Big(Input,varargin{:});

elseif Dimension == 3 && Degree == 3 && MatrixSize== 0 
    [S,Lattice] = computeTriCubicSubdivisionMatrix(Input,varargin{:});

elseif Dimension == 3 && Degree == 3 && MatrixSize== 1 
    [S,Lattice] = computeTriCubicSubdivisionMatrixBig(Input,varargin{:});

end




end