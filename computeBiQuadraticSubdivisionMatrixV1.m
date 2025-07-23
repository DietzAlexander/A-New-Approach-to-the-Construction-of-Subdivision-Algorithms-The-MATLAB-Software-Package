function [S] = computeBiQuadraticSubdivisionMatrixV1(n)
%computeBiQuadraticSubdivisionMatrixV1 compute the subdivision matrix by the 
%standard Doo-Sabin-Algorithm after
%
% Daniel Doo und Malcom Sabin. „Behaviour of Recursive Division Surfaces 
% near Extraordinary Points“. Computer Aided Design 10(6) (1978), S. 356–360.
%
% Input:  n: Amount of points n of the 2-polytope
% Output: S: nxn subdivision matrix
%
%computeBiQuadraticSubdivisionMatrixV1(n) compute the subdivision matrix by the 
%standard Doo-Sabin-Algorithm

%Initializing the output
S=zeros(n,n);

%goes through all rows and colums
for i=1:n
    for j=1:n

        %diagonal entries have a different formula
        if i== j
           S(i,j)=(n+3+2)/(4*n);
        else
            %the distance to the diagonal
            k=j-i;
            if k<0
               k=k+n; 
            end
            S(i,j)=(3+2*cos(2*pi*k/n))/(4*n); 
        end
    end
end

end


