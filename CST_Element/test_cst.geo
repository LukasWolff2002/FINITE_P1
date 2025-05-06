SetFactory("OpenCASCADE");

// Puedo crear 4 puntos 
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {40.0, 0.0, 0.0, 1.0};
Point(3) = {40.0, 10.0, 0.0, 1.0};
Point(4) = {0.0, 10.0, 0.0, 1.0};

// Ahora puedo unir los puntos con lineas
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


// Puedo crear un circulo
Circle(5) = {20, 5, 0, 2.5, 0, 2*Pi};

// Ahora debo crear las superficies,
// Donde su nombre es el que uso para definir la seccion

// Genero un loop a partir de lineas, de esta forma creo el rectangulo
Line Loop(1) = {1,2,3,4};

// Genero un loop del circulo
Curve Loop(2) = {5};

// Ahora defino las superficies

// Superficies del circulo
Plane Surface(1) = {2};

// Superficie del rectangulo - el circulo
Plane Surface(2) = {1,2};

// Ahora defino los grupos fisicos
Physical Surface(1) = {1};  // Superficie central interna
Physical Surface(2) = {2};  // Superficie exterior con hueco interior


// Ahora defino las condiciones de frontera

// En primer lugar defino una linea que sera restringida en X e Y
Physical Line("Restr XY") = {4};  // Línea 4 restringida en X y Y
Physical Line("Fuerza") = {3};    // Línea 2 con fuerza aplicada


