# A-New-Approach-to-the-Construction-of-Subdivision-Algorithms-The-MATLAB-Software-Package
This software package creates subdivision matrices for generalized quadratic and cubic B-spline subdivision. The algorithms produce subdivision matrices for subdivision surfaces as well as subdivision volumes. The subdivision matrices define refinement rules for arbitrary combinatorial structures. Irregular points and edges are supported in all cases. All generated subdivision matrices have a valid eigenstructure: they have a subdominant eigenvalue of 1/2 with multiplicity two (for subdivision surfaces) and multiplicity three (for subdivision volumes). Moreover the central part of the structure forms a convex polytope (or in the cubic case, a set of cubes whose outer faces form a central polytope). The software package also includes a number of additional functions, such as plotting the evaluated B-spline elements or uniformly refining structures of any size.


------------------------------------------------------

To use the software package, it is recommended to take a look at the accompanying

        Manual.pdf

and the two live scripts:

        DemoQuadraticLive.mlx     DemoCubicLive.mlx

------------------------------------------------------


This is the software package accompanying the dissertation 

        Alexander Dietz. Ein neuer Ansatz zur Konstruktion von Subdivisionsalgorithmen
        Dissertation. TU Darmstadt. 2025 DOI: https://doi.org/10.26083/tuprints-00030194

and to the corresponding English translation:

        Alexander Dietz. A New Approach to the Construction of Subdivision Algorithms
        TUprints. Darmstadt. 2025 DOI: https://doi.org/10.26083/tuprints-00030195

Here you will find all the theoretical background.
