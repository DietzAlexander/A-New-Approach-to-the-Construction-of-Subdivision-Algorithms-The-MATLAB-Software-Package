function [S] = computeBiQuadraticSubdivisionMatrixV2(n)
%computeBiQuadraticSubdivisionMatrixV2 compute the subdivision matrix for
%the central 2-polytope after variant 2 of the dissertation.
%
% Input:  n: Amount of points n of the 2-polytope
% Output: S: nxn subdivision matrix
%
%computeBiQuadraticSubdivisionMatrixV2(n) compute the subdivision matrix by 
% the second variant.

%----------------------------------------------------------------------
% Compute the primal points
%----------------------------------------------------------------------

%The edgepoints on the unit circle
Edgepoints=[real(exp(2*pi*1i*((0:n-1))/n))',imag(exp(2*pi*1i*((0:n-1))/n))'];

%The general linear system. Each row is one line orthogonal to the circle
%in and touching the circle in the point edgepoint
GLS=zeros(n,2);

%the right site of the linear syste, 
RS=zeros(n,1);

%compute each line of the linear system and the right site
for i=1:n
    GLS(i,:)=Edgepoints(i,:);
    RS(i)=Edgepoints(i,:)*Edgepoints(i,:)';
end

%The primal points
PrimalPoints=zeros(n,2);

%compute each primal point by solving a linear system consisting out of two
%lines of the GLS
for i=1:n
    if i<n
        PrimalPoints(i,:)=GLS([i,i+1],:)\RS([i,i+1]);
    else
        PrimalPoints(i,:)=GLS([n,1],:)\RS([n,1]);
    end
    
end

%----------------------------------------------------------------------
% Compute the CDV Matrix
%----------------------------------------------------------------------

%compute each entry 
C=zeros(n);
for i=1:n
    for j=1:n
        if i==1
            %if i==1 there are just 2 non diagonal entries unequal to zero
            %the entry n and the entry 2
            if j==n || j==i+1
                C(i,j)=1/sqrt(cross([PrimalPoints(i,:),0],[PrimalPoints(j,:),0])*cross([PrimalPoints(i,:),0],[PrimalPoints(j,:),0])');
            end
        elseif i==n
            %if i==n there are just 2 non diagonal entries unequal to zero
            %the entry n-1 and the entry 1
            if j==i-1 || j==1
                C(i,j)=1/sqrt(cross([PrimalPoints(i,:),0],[PrimalPoints(j,:),0])*cross([PrimalPoints(i,:),0],[PrimalPoints(j,:),0])');
            end
        else
            %for all othe i there are just 2 non diagonal entries unequal to zero
            %the entry i-1 and the entry i+1
            if j==i-1 || j==i+1
                C(i,j)=1/sqrt(cross([PrimalPoints(i,:),0],[PrimalPoints(j,:),0])*cross([PrimalPoints(i,:),0],[PrimalPoints(j,:),0])');
            end
        end
    end
    %conpute the diagon entry
    C(i,i)=-C(i,:)*PrimalPoints/PrimalPoints(i,:);
end

%----------------------------------------------------------------------
% Shifting the spectrum
%----------------------------------------------------------------------


%norming the CDV matrix
C=C./sum(C,2);

%Shifting the spectrum
S=expm((C-eye(n))*log(2));

end
