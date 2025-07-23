PauseTime=4;
clc

Iteration=1;
IterationEnd=26;


fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['This software package is the implementation of the present dissertation:\n\n' ...
    'Alexander Dietz. Ein neuer Ansatz zur Konstruktion von Subdivisionsalgorithmen\n' ...
    'Dissertation. Darmstadt: TU Darmstadt. 2025\n' ...
    'DOI: https://doi.org/10.26083/tuprints-00030194\n\n' ...
    'and the corresponding English translation\n\n' ...
    'Alexander Dietz. A New Approach to the Construction of Subdivision Algorithms\n' ...
    'TUprints. Darmstadt. 2025\n' ...
    'DOI: https://doi.org/10.26083/tuprints-00030195\n\n'])


fprintf('[Press any key to continue.]\n')
pause
clc


%---------------------------------------------------------------
% Quick Use
%---------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['For quick use, there is the function \n\n',...
'[S,Lattice] = computeSubdivisionMatrix(Input,Degree,MatrixSize,varargin) \n\n',...
'This function automatically selects the appropriate algorithm for creating\n', ...
'the subdivision algorithms. The following must be considered:\n\n'])


fprintf(['Input: Must either be a natural number greater than 2\n', ...
    '(for the two-dimensional case) or one of the inputs listed below, \n', ...
    'such as adjacency matrix, edge matrix, or face matrix (for the three-dimensional case).\n\n'])


fprintf("Degree: Must be 2 (or 'quadratic' or 'Quadratic') for the quadratic \n" + ...
    "case or 3 (or 'cubic' or 'Cubic') for the cubic case.\n\n")

fprintf("MatrixSize: Must be 0 (or 'initial' or 'Initial') for initial \n" + ...
    "elements or 1 (or 'big' or 'Big') for larger structures.\n\n");

fprintf("If Variant 1 is to be used in the quadratic case, the command must \n" + ...
    "include the optional input\n\n" + ...
"computeSubdivisionMatrix(Input,Degree,MatrixSize,'mu',value)\n\n"+...
"with 0 <= value < 1/2.\n\n");

fprintf("Also, all optional inputs described below can be used.\n\n" + ...
    "Examples for the command are:\n\n"+...
"S=computeSubdivisionMatrix(5,2,0,'mu',1/4);              2D | quadratic | initial element | V1\n"+...
"S=computeSubdivisionMatrix(6,'quadratic',0);             2D | quadratic | initial element | V2\n"+...
"S=computeSubdivisionMatrix(3,2,'Big','mu',1/4,'Status'); 2D | quadratic | bigger structure | V1\n"+...
"S=computeSubdivisionMatrix(4,2,1,'PreventInputCheck');   2D | quadratic | bigger structure | V2\n"+...
"S=computeSubdivisionMatrix(7,3,'initial');               2D | cubic | initial element \n"+...
"S=computeSubdivisionMatrix(4,'cubic','big');             2D | cubic | bigger structure\n\n")

fprintf("and with\n\n F=\n")
F=[1,2,3,4,5;...
1,2,7,8,6;...
2,3,7,0,0;...
3,4,8,7,0;...
4,5,6,8,0;...
5,1,6,0,0];
disp(F)

fprintf("\nfor the 3D case \n\nS=computeSubdivisionMatrix(F,2,0,'mu',1/4);                                 3D | quadratic | initial element | V1\n"+...
"S=computeSubdivisionMatrix(F,'Quadratic',0);                                3D | quadratic | initial element | V2\n"+...
"[S,Lattice]=computeSubdivisionMatrix(F,2,1,'mu',1/4,'Tolerance',10^(-7));   3D | quadratic | bigger structure | V1\n"+...
"[S,Lattice]=computeSubdivisionMatrix(F,2,'big','MaxIterations',50);         3D | quadratic | bigger structure | V2\n"+...
"[S,Lattice]=computeSubdivisionMatrix(F,3,'Initial');                        3D | cubic | initial element \n"+...
"[S,Lattice]=computeSubdivisionMatrix(F,'cubic','big');                      3D | cubic | bigger structure\n\n")


fprintf('[Press any key to continue.]\n')

pause

clc



%------------------------------------------------------------------------
% Introduction
%------------------------------------------------------------------------
fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf('Welcome to the demo file for the generalized quadratic subdivision algorithm.\n\n')

pause(PauseTime)

fprintf(['The demo will guide you through the key functions, showing their inputs,  \n' ...
    'outputs, and  capabilities. \n \n'])

pause(PauseTime)

fprintf(['After each step, the programm will pause to give you time to review\n' ...
    'the results.\n \n'])



fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Background information
%-------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;
fprintf('This demo provides a brief introduction to the algorithm. \n \n')

pause(PauseTime)

fprintf(['For more details on the code implementation, please refer to \n' ...
    'the corresponding manual.   \n \n'])

pause(PauseTime)

fprintf(['For theoretical background, consult the\n' ...
    'related dissertation:\n\n' ...
    'Alexander Dietz. "Ein neuer Ansatz zur Konstruktion von Subdivisionsalgorithmen"\n' ...
    'Dissertation. TU Darmstadt. 2025 DOI: https://doi.org/10.26083/tuprints-00030194\n \n'])

fprintf(['or the corresponding English translation:\n\n' ...
    'Alexander Dietz. "A New Approach to the Construction of Subdivision Algorithms"\n' ...
    'TUprints. Darmstadt. 2025 DOI: https://doi.org/10.26083/tuprints-00030195\n \n'])

pause(PauseTime)

fprintf(['The code may be used under the\n\n' ...
    'CC-BY 4.0 International\n\n' ...
    'license by citing the corresponding dissertation and its English translation.\n\n'])

pause(PauseTime)

fprintf(['For comments, questions, or to report any errors, you can contact us at any time at:\n\n' ...
    'alexander_dietz@t-online.de \n\n'])







fprintf('[Press any key to continue.]\n')

pause

clc








%------------------------------------------------------------------------
% Two variants
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['For the quadratic case, two variants of the algorithm are available.\n', ...
         'Files for the first variant end with "V1", and files for the second variant end with "V2".\n\n'])

pause(PauseTime)

fprintf("Both variants share a dominant eigenvalue of 1. The subdominant eigenvalue is 1/2 –\n" + ...
        "double for the surface case and triple for the volume case. Both variants can be applied\n" + ...
        "to any combinatorial structure.\n\n")

pause(PauseTime)

fprintf("The first variant has the advantage that the subsubdominant eigenvalue is 1/4.\n\n")

pause(PauseTime)

fprintf("The second variant avoids negative entries and naturally handles semi-irregular\n" + ...
        "volumetric cases through the construction of the subdivision matrix.\n\n")

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Surface Overview
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("Let's start with the surface case.\n\n")

pause(PauseTime)

fprintf(['This case is relatively simple, as both the theory and practical\n' ...
    'applications are well established. \n\n'])

pause(PauseTime)

fprintf(['There are four main files for this case:\n\n' ...
    'computeBiQuadraticSubdivisionMatrixV1(n)\n' ...
    'computeBiQuadraticSubdivisionMatrixV2(n)\n' ...
    'computeBiQuadraticSubdivisionMatrixV1Big(n)\n' ...
    'computeBiQuadraticSubdivisionMatrixV2Big(n)\n\n'])

pause(PauseTime)

fprintf(['Functions without the suffix "Big" compute subdivision  \n' ...
    'matrices for initial elements, refining an \n' ...
    'n-sided polytope into a smaller one.\n\n'])

pause(PauseTime)

fprintf(['Functions with the suffix "Big" compute subdivision \n' ...
    'matrices for a larger area to enable B-spline evaluations,\n' ...
    'refining a ring into a smaller ring of B-spline patches. \n\n'])

pause(PauseTime)

fprintf(['Input: n is a natural number > 2,  indicating the number of\n' ...
    'vertices around the central 2-polytope. \n\n'])

pause(PauseTime)


fprintf('Output: the corresponding subdivision matrix. \n\n')

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% V1 Surface
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("Let's look at an example with n=5. Using the command\n\n" + ... 
    "computeBiQuadraticSubdivisionMatrixV1(5) \n\n" + ...
     "we obtain the subdivision matrix:\n\n")

pause(PauseTime)

n=5;
S1=computeBiQuadraticSubdivisionMatrixV1(n);
disp(S1)

pause(PauseTime)

fprintf('\nThe eigenvalues of this matrix are:\n\n')

sort(eig(S1),'descend')

pause(PauseTime)



fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])


close all;

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));


figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors for an example with n=5 \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));


figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and one refinement \newline  ')

V=S1*V;
plot(V(:,1),V(:,2),'.r','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n') 


B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));


figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and two refinements \newline  ')

V=S1*V;

plot(V(:,1),V(:,2),'.r','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S1*V;

plot(V(:,1),V(:,2),'.r','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% V2 Surface
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("Similar results are obtained for the second variant. Using the command\n\n" + ...
    "computeBiQuadraticSubdivisionMatrixV2(5) \n\n" + ...
    "we obtain the subdivision matrix:\n\n")

pause(PauseTime)

n=5;
S2=computeBiQuadraticSubdivisionMatrixV2(n);
disp(S2)

pause(PauseTime)

fprintf('\nThe eigenvalues of this matrix are: \n\n')

sort(eig(S2),'descend')

pause(PauseTime)
fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

B=S2-1/2*eye(length(S2));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors for an example with n=5 \newline  ')


pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvector results in: \n\n') 


B=S2-1/2*eye(length(S2));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and one refinement \newline  ')

V=S2*V;

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')


B=S2-1/2*eye(length(S2));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and two refinements \newline  ')

V=S2*V;

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S2*V;

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% V1 Surface Big
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['Larger matrices can be computed similarly. Using the command \n\n' ...
    'computeBiQuadraticSubdivisionMatrixV1Big(5) \n\n' ...
    'we obtain the subdivision matrix. Its structure is: \n\n'])

pause(PauseTime)

n=5;
S1Big=computeBiQuadraticSubdivisionMatrixV1Big(n);
close all;

figure
spy(S1Big);
hold on
title('sparse structure of the big subdivision matrix')

fprintf('[Press any key to continue.]\n')
pause


fprintf('\nThe six largest eigenvalues are: \n\n')

EE=sort(eig(S1Big),'descend');
EE(1:6)

pause(PauseTime)
fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors for a big example with n=5 \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and one refinements \newline  ')

V=S1Big*V;

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')

B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and two refinements \newline  ')

V=S1Big*V;

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S1Big*V;

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc



%------------------------------------------------------------------------
% V2 Surface Big
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("Similarly, using the command\n\n" + ...
    "computeBiQuadraticSubdivisionMatrixV2Big(5) \n\n" + ...
     "we obtain the subdivision matrix for variant 2. Its structure is:\n\n")

pause(PauseTime)

n=5;
S2Big=computeBiQuadraticSubdivisionMatrixV2Big(n);
close all;

figure
spy(S2Big);
title('sparse structure of the big subdivision matrix')

fprintf('[Press any key to continue.]\n')
pause


fprintf('\nThe six largest eigenvalues are: \n\n')

EE=sort(eig(S2Big),'descend');
EE(1:6)

pause(PauseTime)
fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

figure
B=S2Big-1/2*eye(length(S2Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors for a big example with n=5 \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

figure
B=S2Big-1/2*eye(length(S2Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and one refinements \newline  ')

V=S2Big*V;
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')

figure
B=S2Big-1/2*eye(length(S2Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and two refinements \newline  ')

V=S2Big*V;
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S2Big*V;
plot(V(:,1),V(:,2),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Overview
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf('We will now begin with the volume case. \n\n')



pause(PauseTime)

fprintf(['There are again 4 main files for this volume case:\n\n' ...
    'computeTriQuadraticSubdivisionMatrixV1(Input,mu,varargin)\n' ...
    'computeTriQuadraticSubdivisionMatrixV2(Input,varargin)\n' ...
    'computeTriQuadraticSubdivisionMatrixV1Big(Input,mu,varargin)\n' ...
    'computeTriQuadraticSubdivisionMatrixV2Big(Input,varargin)\n\n'])

pause(PauseTime)

fprintf(['The functions without the \"Big\" suffix compute subdivision \n' ...
    'matrices for the initial stages of subdivision. These refine a \n' ...
    '3-polytope into a smaller one. \n\n'])

pause(PauseTime)

fprintf(['The functions with the \"Big\" suffix compute subdivision \n' ...
    'matrices for a larger area to support B-spline evaluation.\n' ...
    'Specifically, these functions refine a shell into a smaller \n' ...
    'shell of B-spline cubes. \n\n'])

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Input
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf('We will now look at the general input, which is the parameter \"Input\". \n\n')

pause(PauseTime)

fprintf(['For all four functions, this input represents the combinatorial\n' ...
    'structure of the central polytope.\n\n'])

pause(PauseTime)

fprintf(['There are three possible data formats. The first is the adjacency \n' ...
    'matrix. For the regular case, an example would be: \n\n'])

pause(PauseTime)

A=computePrismFaceMatrix(4);
A=FaceToAdjacencyMatrix(A);

disp(A);

pause(PauseTime)

fprintf(['Alternatively, the input can be a list of faces. \n' ...
    'Each row of this list matrix represents a face, and each entry is a \n' ...
    'vertex of that face. For the regular case, an example would be:\n\n'])

F=computeFaceMatrix(A,false);

disp(F);

pause(PauseTime)

fprintf(['The third option is a list of edges. For the regular case, \n' ...
   'this would look like: \n\n'])

G=graph(A);

disp(G.Edges.EndNodes);

pause(PauseTime)



fprintf(['For the V1 variant, there is an additional input called \"mu\".\n' ...
   'This parameter determines the subsubdominant eigenvalue for the central \n' ...
   'structure and can be chosen by the user. \n\n'])

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Output general
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['Before we look at the optional input parameters, let''s review the output. \n' ...
     'For the functions without the \"Big\" suffix, the output is simply the \n' ...
   'subdivision matrix, which refines a 3-polytope into a smaller one. \n\n'])

pause(PauseTime)

fprintf('For the regular case, we obtain \n\n')

fprintf('A= \n\n')

disp(A);

fprintf('and using the command \n\n')

fprintf(['computeTriQuadraticSubdivisionMatrixV1(A,1/4) \n\n' ...
    'we get: \n\n'])

S=computeTriQuadraticSubdivisionMatrixV1(A,1/4);



disp(S);

pause(PauseTime)

fprintf(['Similarly, for the second variant, using the command \n\n' ...
    'computeTriQuadraticSubdivisionMatrixV2(A) \n\n' ...
    'we get: \n\n'])

S=computeTriQuadraticSubdivisionMatrixV2(A);

pause(PauseTime)

disp(S);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Output V1
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['Next, we will illustrate the algorithms with a more complex example. \n' ...
     'From now on, we will use the face matrix as input.\n\n'])

fprintf('F= \n\n')



F=[1,2,3,4,5;...
1,2,7,8,6;...
2,3,7,0,0;...
3,4,8,7,0;...
4,5,6,8,0;...
5,1,6,0,0];


disp(F)

pause(PauseTime)

fprintf(['For the first variant, using the command\n\n' ...
    'computeTriQuadraticSubdivisionMatrixV1(F,1/4)\n\n' ...
    'we obtain the subdivision matrix: \n\n'])

S1=computeTriQuadraticSubdivisionMatrixV1(F,1/4);


disp(S1);

pause(PauseTime)


fprintf('The eigenvalues of this matrix are: \n\n')

D=sort(eig(S1),'descend');

disp(D);

fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
grid on
title('subdominant eigenvectors for the example with F \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
grid on
title('subdominant eigenvectors and one refinement \newline  ')

V=S1*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
grid on
title('subdominant eigenvectors and two refinements \newline  ')

V=S1*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S1*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Output V2
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['For the second variant, using the command \n\n' ...
    'computeTriQuadraticSubdivisionMatrixV2(F)\n\n' ...
    'we obtain the subdivision matrix: \n\n'])

S1=computeTriQuadraticSubdivisionMatrixV2(F);


disp(S1);

pause(PauseTime)


fprintf('\nThe eigenvalues of this matrix are: \n\n')

D=sort(eig(S1),'descend');

disp(D);

fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
grid on
title('subdominant eigenvectors for the example with F \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
grid on
title('subdominant eigenvectors and one refinement \newline  ')

V=S1*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')

B=S1-1/2*eye(length(S1));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

figure
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
grid on
title('subdominant eigenvectors and two refinements \newline  ')

V=S1*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S1*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc


%------------------------------------------------------------------------
% Volume Output V1 Big
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['The larger matrices are quite similar. Using the command \n\n' ...
    'computeTriQuadraticSubdivisionMatrixV1Big(F,1/4) \n\n' ...
    'we obtain the subdivision matrix with the structure: \n\n'])

pause(PauseTime)

S1Big=computeTriQuadraticSubdivisionMatrixV1Big(F,1/4);
close all;

figure
spy(S1Big);
hold on
title('sparse structure of the subdivision matrix')

fprintf('[Press any key to continue.]\n')
pause


fprintf('\nThe six largest eigenvalues of this matrix are: \n\n')

EE=sort(eig(S1Big),'descend');
EE(1:6)

pause(PauseTime)
fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

figure
B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors for a big example with F \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

figure
B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and one refinement \newline  ')

V=S1Big*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')

figure
B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and two refinements \newline  ')

V=S1Big*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S1Big*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Output V2 Big
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['Using the command \n\n' ...
    'computeTriQuadraticSubdivisionMatrixV2Big(F) \n\n' ...
   'we obtain the subdivision matrix for variant 2. The structure of the\n' ...
    'matrix is: \n\n'])

pause(PauseTime)

S1Big=computeTriQuadraticSubdivisionMatrixV2Big(F);
close all;

figure
spy(S1Big);
title('sparse structure of the subdivision matrix')

fprintf('[Press any key to continue.]\n')
pause


fprintf('\nThe six largest eigenvalues of this matrix are: \n\n')

EE=sort(eig(S1Big),'descend');
EE(1:6)

pause(PauseTime)
fprintf(['\n' ...
    'The subdominant eigenvectors can be plotted as follows: \n\n'])

close all;

figure
B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors for a big example with F \newline  ')

pause(PauseTime)

fprintf('Multiplying the subdivision matrix by the eigenvectors results in: \n\n')

figure
B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and one refinement \newline  ')

V=S1Big*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

pause(PauseTime)

fprintf('Once more: \n\n')

figure
B=S1Big-1/2*eye(length(S1Big));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));

plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0 0.4470 0.7410]);
hold on
axis equal
axis off
title('subdominant eigenvectors and two refinements \newline  ')

V=S1Big*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.8500 0.3250 0.0980]);

V=S1Big*V;
plot3(V(:,1),V(:,2),V(:,3),'.','Markersize',25,'Color',[0.9290 0.6940 0.1250]);

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Volume Output Lattice
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['For the functions with the \"Big\" suffix, there is an additional output parameter\n' ...
    'called \"Lattice\". This parameter encodes the structure of the control points \n' ...
    'or, more precisely, the underlying structure. Let n represent the number of control \n' ...
     'points in the structure. Each column of \"Lattice\" then has size n^2, since each \n' ...
     'column represents an adjacency matrix of a polytope in this structure. \n' ...
     'The original adjacency matrices can be reconstructed using the command \n\n' ...
    'reshape(Lattice(:,i),[n,n]) \n\n'])

pause(PauseTime)

fprintf(['The structure of the \"Lattice\" matrix for our example is \n' ...
    '(applies to both variants): \n\n'])


[S1Big,Lattice]=computeTriQuadraticSubdivisionMatrixV2Big(F);

close all

figure
spy(Lattice)
title('sparse structure of lattice matrix \newline ')

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Additional Commands Status
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf(['Now we will take a look at the additional input parameters.\n' ...
    'These parameters are the same for all four variants. Therefore, \n' ...
   'we will focus on the output for variant 1. \n\n'])

pause(PauseTime)

fprintf("The first additional input argument is 'Status'. If this \n" + ...
    "parameter is included, the algorithm will display the status \n"+ ...
    "of the construction process. Using the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Status') \n\n" + ...
    "we get:\n\n")

fprintf('[Press any key to continue.]\n')
pause
clc

computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Status');

fprintf('\n [Press any key to continue.]\n')
pause
clc

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("Similarly, for the big function, using the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1Big(F,1/4,'Status') \n\n" + ...
    "we get:\n\n")

fprintf('[Press any key to continue.]\n')
pause
clc

computeTriQuadraticSubdivisionMatrixV1Big(F,1/4,'Status');

fprintf('\n [Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Additional Commands Status and tolerance
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("The second additional input argument is 'PreventInputCheck'.\n" + ...
    "If this parameter is included, the algorithm will skip the check \n"+ ...
    "to ensure that the input represents a planar 3-connected graph. \n\n" + ...
    "The command is\n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'PreventInputCheck')\n\n")

pause(PauseTime)

fprintf("The third additional input argument is 'Tolerance'. \n" + ...
    "This parameter allows you to manually set the tolerance for the \n"+ ...
    "numerical parts of the algorithm. The default value is 10^(-13).\n" + ...
    "The command is \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Tolerance',10^(-8)) \n\n")

pause(PauseTime)

fprintf("You can also manually set the maximum number of iterations\n" + ...
    "for the numerical parts using the 'MaxIterations' parameter. \n" + ...
    "The default value is 100. To change it, use the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'MaxIterations',50) \n\n")

pause(PauseTime)

fprintf('Of course, these options also apply to the functions with the "Big" suffix.\n\n')

fprintf('\n[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Additional Commands Visualization
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("Another input parameter is 'Visualization'.\n" + ...
    "With this parameter, various graphics are generated during the \n" + ...
    "construction process. Use the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Visualization') \n\n")

fprintf('\n[Press any key to continue.]\n')
pause

close all
computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Visualization');

fprintf("\nA total of 13 figures are created. The background and explanations \n" + ...
    "for these figures can be found in the corresponding dissertation. \n\n")

pause(PauseTime)

fprintf(['For the functions with the "Big" suffix, only one figure is created. It can be\n', ...
         'generated with the command\n\n', ...
         'computeTriQuadraticSubdivisionMatrixV1Big(F,1/4,''Visualization'')\n\n'])

fprintf('\n[Press any key to continue.]\n')
pause

figure
computeTriQuadraticSubdivisionMatrixV1Big(F,1/4,'Visualization');

fprintf('\n[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Additional Commands Visualization - 2
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("These figures can be manipulated in different ways. For example,\n" + ...
    "with the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Visualization','DrawAfterXKites',4) \n\n" + ...
    "the initial iterative figure will display faster as four kites are drawn \n" + ...
    "in one step.\n\n")

pause(PauseTime)

fprintf("Using the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Visualization','View',[26,5]) \n\n" + ...
    "the viewing angle of the figures can be adjusted. \n\n")

pause(PauseTime)

fprintf("And with the command \n\n" + ...
    "computeTriQuadraticSubdivisionMatrixV1(F,1/4,'Visualization','PrintImages','Test','.png') \n\n" + ...
    "the graphics are exported to .png files with the prefix 'Test' and the suffix '.png'. \n\n")

fprintf('\n[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Quality of life - plot
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("There are also two quality-of-life functions. The first one is a\n" + ...
    "plot function. We start with an example using the following commands\n\n" + ...
    "[S,Lattice]=computeTriQuadraticSubdivisionMatrixV1Big(F,1/4);\n" + ...
    "B=S-1/2*eye(length(S));\n" + ...
    "B=B'*B;\n" + ...
    "[V,D]=eig(B);\n" + ...
    "D=diag(D);\n" + ...
    "a=abs(D)<10^(-10);\n" + ...
    "V=V(:,a);\n" + ...
    "V=real(orth(V));\n" + ...
    "Controlpoints=V; \n\n")

pause(PauseTime)

fprintf("With this, we have a lattice and a set of control points \n" + ...
    "belonging to the eigenshell of the structure. The command\n\n" + ...
    "plotTriQuadraticBSplineLattice(Controlpoints,Lattice) \n\n" + ...
    "plots all B-spline cubes that can be evaluated within this structure.\n" + ...
    "The output is:\n")

fprintf('\n[Press any key to continue.]\n')
pause


[S,Lattice]=computeTriQuadraticSubdivisionMatrixV1Big(F,1/4);
B=S-1/2*eye(length(S));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));
Controlpoints=V;

close all
figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice)
title('plot of the evaluable B-spline cubes')

fprintf('\n[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Quality of life - plot - 2
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("The plot command has two optional input parameters. The first \n" + ...
    "one is 'PointsPerCube'. This parameter allows you to specify how \n" + ...
    "many points each cube should be evaluated in each direction. Using \n" + ...
    "the command \n\n" + ...
    "plotTriQuadraticBSplineLattice(Controlpoints,Lattice,'PointsPerCube',2)\n\n" + ...
    "we get:")

fprintf('\n[Press any key to continue.]\n')
pause

figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice,'PointsPerCube',2)
title('plot with only 2x2x2 points per cube')

fprintf('\n[Press any key to continue.]\n')
pause

fprintf("\nThe second optional parameter is 'ColorTable', which allows \n" + ...
    "you to adjust the colors of the cubes. For example, we define:\n\n" + ...
    "CT=\n")

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

CT=[ColorTable(:,[2,1,3])];

disp(CT)

fprintf("\nUsing the command \n\n" + ...
   "plotTriQuadraticBSplineLattice(Controlpoints,Lattice,'ColorTable',CT)\n\n" + ...
   "we get:")


fprintf('\n[Press any key to continue.]\n')
pause

figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice,'ColorTable',CT)
title('plot with different colors')

fprintf('[Press any key to continue.]\n')
pause
clc

%------------------------------------------------------------------------
% Quality of life uniform refinement
%------------------------------------------------------------------------

fprintf("Part "+Iteration +" of "+IterationEnd +"\n\n")
Iteration=Iteration+1;

fprintf("The second quality-of-life function is  uniform refinement. \n" + ...
    "This function takes a matrix of control points and a lattice matrix \n" + ...
    "and refines the entire structure uniformly. The available functions are: \n\n" + ...
    "RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,mu)\n" + ...
    "RefineTriQuadraticLatticeUniformV2(Controlpoints,Lattice) \n\n")

pause(PauseTime)

fprintf("The output of these functions is a refined matrix of control points \n" + ...
    "and the updated lattice matrix. Additional output options are available \n" + ...
    "for extracting volumes out of volumes, faces, edges, and points from the lattice (as columns). \n" + ...
    "The functions can also return which vertices are inner vertices and their valences.\n\n")

pause(PauseTime)

fprintf("A simple example is as follows:\n\n" + ...
    "[S,Lattice]=computeTriQuadraticSubdivisionMatrixV1Big(F,1/4);\n" + ...
    "B=S-1/2*eye(length(S));\n" + ...
    "B=B'*B;\n" + ...
    "[V,D]=eig(B);\n" + ...
    "D=diag(D);\n" + ...
    "a=abs(D)<10^(-10);\n" + ...
    "V=V(:,a);\n" + ...
    "V=real(orth(V));\n" + ...
    "Controlpoints=V; \n"+ ...
    "plotTriQuadraticBSplineLattice(Controlpoints,Lattice)\n" + ...
    "[Controlpoints,Lattice]=RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,1/4)\n" + ...
    "figure\n" + ...
    "plotTriQuadraticBSplineLattice(Controlpoints,Lattice) \n\n" + ...
    "This produces: \n\n")

fprintf('[Press any key to continue.]\n')
pause

close all

[S,Lattice]=computeTriQuadraticSubdivisionMatrixV1Big(F,1/4);
B=S-1/2*eye(length(S));
B=B'*B;
[V,D]=eig(B);
D=diag(D);
a=abs(D)<10^(-10);
V=V(:,a);
V=real(orth(V));
Controlpoints=V;
figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice);
title('start subdivision volume')

[Controlpoints,Lattice]=RefineTriQuadraticLatticeUniformV1(Controlpoints,Lattice,1/4);
figure
plotTriQuadraticBSplineLattice(Controlpoints,Lattice);
title('subdivision volume after refinement')

fprintf('\n[Press any key to continue.]\n')
pause

fprintf("\nThis function can also handle lattices containing more than one irregular vertex.\n" + ...
    "Here is an example: \n\n")

fprintf('[Press any key to continue.]\n')
pause

close all
ExampleQuadraticUniformRefinement

fprintf('[Press any key to continue.]\n')
pause
clc

fprintf("Thanks for watching, this demo ends now.\n")