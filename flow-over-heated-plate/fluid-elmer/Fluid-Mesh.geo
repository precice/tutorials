//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 0.5, 0, 1.0};
//+
Point(4) = {0, 0.5, 0, 1.0};
//+
Point(5) = {-0.5, 0, 0, 1.0};
//+
Point(6) = {-0.5, 0.5, 0, 1.0};
//+
Point(7) = {3, 0, 0, 1.0};
//+
Point(8) = {3, 0.5, 0, 1.0};
//+
Line(1) = {5, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 7};
//+
Line(4) = {7, 8};
//+
Line(5) = {8, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 6};
//+
Line(8) = {6, 5};
//+
Line Loop(1) = {6, 7, 8, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("fluid") = {1};
//+
Physical Line("Inlet") = {8};
//+
Physical Line("Outlet") = {4};
//+
Physical Line("Coupling_Interface") = {2};
//+
Physical Line("Pipe_Boundary") = {1, 7, 6, 5, 3};
//+
Transfinite Surface {1} = {5, 7, 8, 6};
//+
Transfinite Line {8, 4, 1, 7} = 5 Using Progression 1;
//+
Transfinite Line {2, 6} = 10 Using Progression 1;
//+
Transfinite Line {3, 5} = 20 Using Progression 1;
