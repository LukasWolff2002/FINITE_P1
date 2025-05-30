SetFactory("OpenCASCADE");

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {2800.0, 0.0, 0.0, 1.0};
Point(3) = {2800.0, 600.0, 0.0, 1.0};
Point(4) = {0.0, 600, 0.0, 1.0};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Surface('Concrete') = {1}; 

Physical Line("Restr XY") = {4};  // Línea 4 restringida en X y Y
Physical Line("Fuerza 1") = {2};  // Línea 4 restringida en X y Y


