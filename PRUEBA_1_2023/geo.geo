SetFactory("OpenCASCADE");

Point(1) = {1000.0, 0.0, 0.0, 1.0};
Point(2) = {1600.0, 0.0, 0.0, 1.0};
Point(3) = {1600.0, 1200.0, 0.0, 1.0};
Point(4) = {2600.0, 1200.0, 0.0, 1.0};
Point(5) = {2600.0, 1600.0, 0.0, 1.0};
Point(6) = {0.0, 1600.0, 0.0, 1.0};
Point(7) = {0.0, 1200.0, 0.0, 1.0};
Point(8) = {1000.0, 1200.0, 0.0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};

Line Loop(1) = {1,2,3,4,5,6,7,8};

Plane Surface(1) = {1};

Physical Surface(1) = {1};  // Superficie central interna

Physical Line("Restr 3") = {1};  // Línea 4 restringida en X y Y
Physical Line("Restr 2") = {4};  // Línea 4 restringida en X y Y
Physical Line("Restr 1") = {6};  // Línea 4 restringida en X y Y
Physical Line("Force") = {5};  // Línea 4 restringida en X y Y
