function globefun3(R,lat,lon,bool_point,lc)
%GLOBEFUN3 plots a single point on the sphere with a reference globe
%
% INPUT
%   R           radius of sphere
%   lat         latitude (deg)  (theta = 90 - lat)
%   lon         longitude (deg) (lon = phi)
%   bool_point  plot the specified point
%   lc          line style/color
%
% EXAMPLE: figure; globefun3(1,60,40,1,'r');
%
% see globe.m, globefunN.m
% called by test_globefun.m
%

deg = 180/pi;
lw = 1.0;
hold on;

th = (90-lat)/deg;
ph = lon/deg;

% (1) xy-plane circle (equator)
% (2) yz-plane
% (3) xz-plane (prime meridian)
gpts = 100;
tt = linspace(0, 2*pi, gpts);
plot3(R*cos(tt), R*sin(tt), zeros(1,gpts), lc, 'linewidth', lw);
plot3(zeros(1,gpts), R*cos(tt), R*sin(tt), lc, 'linewidth', lw);
plot3(R*cos(tt), zeros(1,gpts), R*sin(tt), lc, 'linewidth', lw);

% axes
fac = 1.15; lw2 = 0.5;
plot3(fac*[-R R],[0 0],[0 0],'k', 'linewidth', lw2);
plot3([0 0],fac*[-R R],[0 0],'k', 'linewidth', lw2);
plot3([0 0],[0 0],fac*[-R R],'k', 'linewidth', lw2);

if bool_point == 1
	% plot the specified point P(r,th,phi)
	msize = 20; lw3 = lw+0.5;
    % [X,Y,Z] = SPH2CART(az,el,r)
    [X,Y,Z] = sph2cart(ph,pi/2-th,R);
    [X2,Y2,Z2] = sph2cart(ph,0,R);

	plot3(X,Y,Z,[lc '.'],'markersize',msize);
	plot3([0 X],[0 Y],[0 Z],lc, 'linewidth', lw3);
	
    leps = 0.1;
    if and( lat >= -90+leps, lat <= 90-leps)
        xzplane = [R*cos(tt) ; zeros(1,gpts) ; R*sin(tt)];
        rotmat = [cos(ph) -sin(ph) 0 ; sin(ph) cos(ph) 0; 0 0 1];
        rotxz = rotmat*xzplane;
        plot3(rotxz(1,:), rotxz(2,:), rotxz(3,:),lc,'linewidth',lw);
        aa = lat/deg;
        plot3(R*cos(aa)*cos(tt), R*cos(aa)*sin(tt), R*sin(aa)*ones(gpts, 1),lc,'linewidth', lw);
        
        % plots some extra lines for reference
        plot3([0 X2],[0 Y2],[0 0],lc, 'linewidth', lw2);
        plot3([X X],[Y Y],[0 Z],lc, 'linewidth', lw2);
        plot3([0 X],[0 Y],[Z Z],lc, 'linewidth', lw2);
    end
    %title(sprintf('\\phi = %.2f, \\theta = %.2f',ph*deg,th*deg));
end

% origin
plot3(0,0,0,'kx');
view(100, 10); axis equal; box on;
xlabel('x'); ylabel('y'); zlabel('z');

%==========================================================================
