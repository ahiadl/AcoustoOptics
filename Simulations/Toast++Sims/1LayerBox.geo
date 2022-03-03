res = 0.5;
yWidth = 30;
zWidth = 30;
xMid = 20;
xEnd = 20;

Point(1) = {0, -yWidth,  zWidth, res};
Point(2) = {0,  yWidth,  zWidth, res};
Point(3) = {0,  yWidth, -zWidth, res};
Point(4) = {0, -yWidth, -zWidth, res};
Point(5) = {xMid, -yWidth,  zWidth, res};
Point(6) = {xMid,  yWidth,  zWidth, res};
Point(7) = {xMid,  yWidth, -zWidth, res};
Point(8) = {xMid, -yWidth, -zWidth, res};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9)  = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};


Line Loop(13) = {1, 2, 3, 4};
Plane Surface(13) = {13};
Line Loop(14) = {9, 10, 11, 12};
Plane Surface(14) = {14};
Line Loop(15) = {1, 6, -9, -5};
Plane Surface(15) = {15};
Line Loop(16) = {6, 10, -7, -2};
Plane Surface(16) = {16};
Line Loop(17) = {-3, 7, 11, -8};
Plane Surface(17) = {17};
Line Loop(18) = {-5, -4, 8, 12};
Plane Surface(18) = {18};

Surface Loop(19) = {13, 15, 16, 17, 18, 14};
Volume(19) = {19};