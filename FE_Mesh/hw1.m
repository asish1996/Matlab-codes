clear all
clc

Lx = 1.5;
Ly = 1;
nx = 4;
ny = 3;

L = Lx * nx;
B = Ly * ny;

%corners = [-L/2,-B/2;L/2,-B/2;L/2,B/2;-L/2,B/2];
corners = [0,0;5,0;5,6;0,5];
%% Defining boundaries
b1x = [];
b1y = [];
for i = 1:(nx+1)
    b1x(i) = corners(1,1) + (i-1)*((corners(2,1)-corners(1,1))/nx);
    b1y(i) = corners(1,2) + (i-1)*((corners(2,2)-corners(1,2))/nx);
end
b1 = [b1x',b1y'];    

b2x = [];
b2y = [];
for i = 1:(ny+1)
     b2x(i) = corners(2,1) + (i-1)*((corners(3,1)-corners(2,1))/ny);
     b2y(i) = corners(2,2) + (i-1)*((corners(3,2)-corners(2,2))/ny);
end
b2 = [b2x',b2y'];

b3x = [];
b3y = [];
for i = 1:(nx+1)
    b3x(i) = corners(3,1) + (i-1)*((corners(4,1)-corners(3,1))/nx);
    b3y(i) = corners(3,2) + (i-1)*((corners(4,2)-corners(3,2))/nx);
end
b3 = [b3x',b3y'];

b4x = [];
b4y = [];
for i = 1:(ny+1)
    b4x(i) = corners(4,1) + (i-1)*((corners(1,1)-corners(4,1))/ny);
    b4y(i) = corners(4,2) + (i-1)*((corners(1,2)-corners(4,2))/ny);
end
b4 = [b4x',b4y'];

%% Node coordinates
Xcoords = zeros(ny+1,nx+1);

for i = 1:(ny+1)
    for j = 1:(nx+1)
        Xcoords(i,j) = b4(ny+2-i,1) + (j-1)*((b2(i,1) - b4(ny+2-i,1))/nx);
    end
end
Xcoor = Xcoords';
Xcoord_vect = Xcoor(:);

Ycoords = zeros(nx+1,ny+1);
for i = 1:(nx+1)
    for j = 1:(ny+1)
        Ycoords(i,j) = b1(i,2) + (j-1)*((b3(nx+2-i,2) - b1(i,2))/ny);
    end
end
Ycoord_vect = Ycoords(:);
Node_coords = [Xcoord_vect,Ycoord_vect];

%% Connectivity Matrix
connec1 = [1:((nx+1)*(ny+1)-(nx+2))]';
connec2 = [2:((nx+1)*(ny+1)-(nx+1))]';
connec3 = [nx+3:(nx+1)*(ny+1)]';
connec4 = [nx+2:((nx+1)*(ny+1)-1)]';
for i = 1:ny-1
    connec1(5*i - (i-1)) = [];
    connec2(5*i - (i-1)) = [];
    connec3(5*i - (i-1)) = [];
    connec4(5*i - (i-1)) = [];
end
connec = [connec1,connec2,connec3,connec4];
%% Boundary matrix
e11 = [1:nx]';
e12 = [2:nx+1]';
e21 = [nx+1:nx+1:(nx+1)*(ny+1) - nx]';
e22 = [2*(nx+1):nx+1:(nx+1)*(ny+1)]';
e31 = [(nx+1)*(ny+1):-1:(nx+1)*(ny+1)-(nx-1)]';
e32 = [(nx+1)*(ny+1)-1:-1:(nx+1)*(ny+1)-nx]';
e41 = [(nx+1)*(ny+1)-nx:-nx-1:nx+2]';
e42 = [(nx+1)*(ny+1)-2*nx-1:-nx-1:1]';
edge1 = ones(nx,1);
edge2 = 2*ones(ny,1);
edge3 = 3*ones(nx,1);
edge4 = 4*ones(ny,1);
bound = [e11,e12,edge1;e21,e22,edge2;e31,e32,edge3;e41,e42,edge4];

%% Plotting mesh

for i = 1:nx+1
    plot ([b1(i,1) b3(nx+2-i,1)],[b1(i,2) b3(nx+2-i,2)])
    hold on
end

for i = 1:ny+1
    plot ([b4(ny+2-i,1) b2(i,1)],[b4(ny+2-i,2) b2(i,2)])
    hold on
end

dlmwrite('bound.txt',bound)

%Please write whichever ever data you would like to extract to .txt file in
%dlmwrite function
%Node_coords array gives us location of nodes
%connec array gives the connectivity
%bound arrary gives the boundaries in specified format of the question