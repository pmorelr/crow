// dimensions

x1=-1.0;  	// first (inner) rectangle			
x2=0.;
y1=1.0;
y2=-1.0;					


X1=-2.;  	// second (box boundary) rectangle			
X2=0.05;
Y1=2.;
Y2=-2.;


xc=-0.25;
yc=0.0;
rc=0.01;

Xc=0.25;
Yc=0.0;
Rc=0.01;


// resolutions
//resolutionc=0.025;
r1=0.02;
r2=0.04;
r3=0.2;


// first (inner) rectangle
Point(1) = {x1,y1,0,r1};
Point(2) = {x2,y1,0,r1};
Point(3) = {x2,y2,0,r1};
Point(4) = {x1,y2,0,r1};
// second rectangle (box boundary)
Point(5) = {X1,Y1,0,r3};
Point(6) = {X2,Y1,0,r3};
Point(7) = {X2,y1,0,r1};
Point(8) = {X2,y2,0,r1};
Point(9) = {X2,Y2,0,r3};
Point(10) = {X1,Y2,0,r3};

// first cylinder
//Point(9) = {xc,yc,0,r1};
//Point(10) = {xc-rc,yc,0,r1};
//Point(11) = {xc,yc+rc,0,r1};
//Point(12) = {xc+rc,yc,0,r1};
//Point(13) = {xc,yc-rc,0,r1};
// second cylinder
//Point(14) = {Xc,Yc,0,r1};
//Point(15) = {Xc-Rc,Yc,0,r1};
//Point(16) = {Xc,Yc+Rc,0,r1};
//Point(17) = {Xc+Rc,Yc,0,r1};
//Point(18) = {Xc,Yc-Rc,0,r1};


// first (inner) rectangle
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// second rectangle (box boundary)
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10)= {10,5};

// first cylinder
//Circle(9) = {10,9,11};
//Circle(10) = {11,9,12};
//Circle(11) = {12,9,13};
//Circle(12) = {13,9,10};
// second cylinder
//Circle(13) = {15,14,16};
//Circle(14) = {16,14,17};
//Circle(15) = {17,14,18};
//Circle(16) = {18,14,15};


Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8,9,10};
//Line Loop(3) = {9,10,11,12};
//Line Loop(4) = {13,14,15,16};


Plane Surface(1) = {1};	// inner rectangle
Plane Surface(2) = {2,1};	// second rectangle minus inner rectangle
//Plane Surface(1) = {1,3,4};	// inner rectangle
//Plane Surface(2) = {1,2,3,4};	// second rectangle minus inner rectangle

Physical Surface(1) = {1,2};  // add all to form one physical domain
Physical Line(1) = {1};	// top
Physical Line(2) = {2};	// right
Physical Line(3) = {3};	// bottom
Physical Line(4) = {4};	// left


