% creat great circle cross section with mesh grids in geocentric coordinate

clear all;close all;clc
%% parameters

% referece spheroid
E = referenceEllipsoid('wgs84');
ecc = E.Eccentricity;

% Normalization dimesions used in SPECFEM
R_earth = 6371; % km

% begin/end points on the earch surface (geodetic degree)
%YN.GYA
lon0 = 106.664;
lat0 = 26.459;
%HYPO
lon1 = 139.23;
lat1 = 46.99;
r0 = R_earth;
r1 = R_earth-900;

% mesh size
na = 100;
nr = 100;

% output
outfile = 'section.xyz';

%% useful formula

% degree to radian 
deg2rad = pi/180;

% rotation matrix: rotate around unit vector (n) by angle (theta)
R = @(n,theta) cos(theta)*eye(3) ...
             + (1-cos(theta))*reshape(n,[3,1])*reshape(n,[1 3]) ...
             + sin(theta)*[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
           
% unit vector specified by (theta,phi)
V = @(theta,phi) [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

%% get grid steps on the radial and radius direction

% get geocentric (theta,phi)
theta0 = pi/2-geodetic2geocentricLat(ecc,lat0*deg2rad);
theta1 = pi/2-geodetic2geocentricLat(ecc,lat1*deg2rad);

phi0 = lon0*deg2rad;
phi1 = lon1*deg2rad;

% get the normal of the cross-section
v0 = V(theta0,phi0);
v1 = V(theta1,phi1);
n = cross(v0,v1); 
n = n/norm(n);
           
% the angle between v0 and v1
a01 = acos(dot(v0,v1));

% get the radial angle steps 
alpha = linspace(0,a01,na);
va = zeros(3,na);
for i = 1:na
  va(:,i) = R(n,alpha(i))*v0;
end

% get the radius steps
r = linspace(r0,r1,nr);
r = reshape(r,[nr,1]);

%% form mesh grid

xx = r*va(1,:)/R_earth;
yy = r*va(2,:)/R_earth;
zz = r*va(3,:)/R_earth;
rr = r(:,ones(na,1));
aa = alpha(ones(nr,1),:);

%% save mesh grids

save('mesh.mat');

%% output

fid = fopen(outfile,'w');
for i = 1:numel(xx)
  fprintf(fid,'%10.8f %10.8f %10.8f\n',xx(i),yy(i),zz(i));
end
fclose(fid);

%% END