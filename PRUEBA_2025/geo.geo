SetFactory("OpenCASCADE");

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {2800.0, 600.0, 0.0, 1.0};
Point(3) = {2800.0, 1600.0, 0.0, 1.0};
Point(4) = {610.0, 1600, 0.0, 1.0};
Point(5) = {605.0, 1600-200, 0.0, 1.0};
Point(6) = {600.0, 1600, 0.0, 1.0};
Point(7) = {0.0, 1600, 0.0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,1};

Line Loop(1) = {1,2,3,4,5,6,7};

Plane Surface(1) = {1};

Physical Surface('Concrete') = {1}; 

Physical Line("Restr XY") = {7};  // Línea 4 restringida en X y Y
Physical Line("Fuerza 1") = {6};  // Línea 4 restringida en X y Y
Physical Line("Fuerza 2") = {3};  // Línea 4 restringida en X y Y

