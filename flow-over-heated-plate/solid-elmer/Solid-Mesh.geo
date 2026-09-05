Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, -0.25, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {1, -0.25, 0, 1.0};
//+
Line(1) = {2, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 2};
//+
Line Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Line("Plate_Bottom") = {1};
//+
Physical Line("Plate_Sides") = {4, 2};
//+
Physical Line("Coupling_Interface") = {3};
//+
Transfinite Surface {1} = {2, 4, 3, 1};
//+
Transfinite Line {4, 2} = 10 Using Progression 1;
//+
Transfinite Line {1, 3} = 10 Using Progression 1;
//+
Physical Surface("Plate") = {1};
