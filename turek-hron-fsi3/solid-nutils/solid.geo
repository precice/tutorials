// Geometry file for the turek.py example.
//
// This is a generalized description of the setup defined by Turek and Hron,
// which requires the following lenghts to be supplied externally, using the
// numbers argument of mesh.gmsh or the -setnumber switch when invoking the
// gmsh application directly:
//
// - x_center: horizontal position of the cylinder measured from the left edge
// - y_center: vertical position of the cylinder measured from the bottom edge
// - cylinder_radius: radius of the cylinder
// - structure_length: length of the elastic structure measured from the cylinder wall
// - structure_thickness: thickness of the elastic structure
// - elemsize: mesh element size in the solid
//
// The parameterization matches largely that of Table 1 of Turek and Hron 2006,
// with the main difference that reference point A cannot be independently
// placed but are always located at the tip of the elastic structure and the
// leading edge of the cylinder, respectively.

SetFactory("OpenCASCADE");

Rectangle(1) = {x_center, y_center - structure_thickness/2, 0, cylinder_radius + structure_length, structure_thickness, 0};
Disk(2) = {x_center, y_center, 0, cylinder_radius};
BooleanIntersection(3) = { Surface{1}; }{ Surface{2}; };
BooleanDifference(4) = { Surface{1}; }{ Surface{2}; };
A = newp; Point(A) = {x_center + cylinder_radius + structure_length, y_center, 0};

// At this point surface 3 (rectangle in cylinder) and 4 (solid domain) are
// non-overlapping. Gmsh promises that the boolean fragments operation with
// deletion will reuse the surface IDs for the new objects.

_() = BooleanFragments{ Surface{3,4}; Point{A}; Delete; }{};

// Fragments deduplicates boundary segments, which means that we can now
// perform boolean operations on the index sets.

tmp() = Abs(Boundary{ Surface{3}; });
bnd_structure() = Abs(Boundary{ Surface{4}; });
bnd_cylinder = bnd_structure();
bnd_structure -= tmp();
bnd_cylinder -= bnd_structure();

Physical Surface("all") = {4};
Physical Line("fluid") = {bnd_structure()};
Physical Line("cylinder") = {bnd_cylinder()};
Physical Point("A") = {A};

// The element size is set equal to elemsize throughout the domain.

Mesh.MeshSizeMin = elemsize;
Mesh.MeshSizeMax = elemsize;
