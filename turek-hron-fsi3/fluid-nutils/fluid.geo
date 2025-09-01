// Geometry file for the turek.py example.
//
// This is a generalized description of the setup defined by Turek and Hron,
// which requires the following lenghts to be supplied externally, using the
// numbers argument of mesh.gmsh or the -setnumber switch when invoking the
// gmsh application directly:
//
// - channel_length: length of the fluid domain
// - channel_height: height of the fluid domain
// - x_center: horizontal position of the cylinder measured from the left edge
// - y_center: vertical position of the cylinder measured from the bottom edge
// - cylinder_radius: radius of the cylinder
// - structure_length: length of the elastic structure measured from the cylinder wall
// - structure_thickness: thickness of the elastic structure
// - min_elemsize: mesh element size at the solid/fluid interface
// - max_elemsize: mesh element size at the channel wall and far field
//
// The parameterization matches largely that of Table 1 of Turek and Hron 2006,
// with the main difference that reference points A and B cannot be
// independently placed but are always located at the tip of the elastic
// structure and the leading edge of the cylinder, respectively.

SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, channel_length, channel_height};
Rectangle(2) = {x_center, y_center - structure_thickness/2, 0, cylinder_radius + structure_length, structure_thickness, 0};
Disk(3) = {x_center, y_center, 0, cylinder_radius};
BooleanDifference(4) = { Surface{2}; }{ Surface{3}; };
BooleanDifference(5) = { Surface{1}; }{ Surface{2,3}; };
A = newp; Point(A) = {x_center + cylinder_radius + structure_length, y_center, 0};
B = newp; Point(B) = {x_center - cylinder_radius, y_center, 0};

// At this point surface 3 (cylinder), 4 (solid domain) and 5 (fluid domain) are
// non-overlapping. Gmsh promises that the boolean fragments operation with
// deletion will reuse the surface IDs for the new objects.

_() = BooleanFragments{ Surface{3,4,5}; Point{A,B}; Delete; }{};

// Fragments deduplicates boundary segments, which means that we can now
// perform boolean operations on the index sets.

bnd_cylinder() = Abs(Boundary{ Surface{3}; });
bnd_structure() = Abs(Boundary{ Surface{4}; });
tmp = bnd_structure();
bnd_structure -= bnd_cylinder();
bnd_cylinder -= tmp();
bnd_fluid() = Abs(Boundary{ Surface{5}; });
bnd_fluid -= bnd_structure();
bnd_fluid -= bnd_cylinder();

// After subtracting the inner boundaries, only the four boundary segments of
// rectangle 1 remain in bnd_fluid, and we are going to assume that they are
// ordered bottom, right, top, left.

Physical Surface("fluid") = {5};
Physical Line("inlet") = {bnd_fluid(3)};
Physical Line("outlet") = {bnd_fluid(1)};
Physical Line("wall") = {bnd_fluid(0), bnd_fluid(2)};
Physical Line("cylinder") = {bnd_cylinder()};
Physical Line("structure") = {bnd_structure()};
Physical Point("A") = {A};
Physical Point("B") = {B};

// The element size is set to be uniformly min_elemsize inside the elastic
// structure, and grow linearly to max_elemsize in the fluid domain over a
// distance of half the channel height.

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Field[1] = Distance;
Field[1].SurfacesList = {3,4};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMin = 0;
Field[2].DistMax = channel_height/2;
Field[2].SizeMin = min_elemsize;
Field[2].SizeMax = max_elemsize;
Background Field = 2;
