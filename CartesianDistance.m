function [CartDistance] = CartesianDistance(Point1, Point2)
% compute Cartesian distnse in [km] between 2 points on spherical Earth
% input angles Lat, Long in degrees
% by Alexandr Sokolov
Re = 6378100; % Radius of Earth, [m]

%                                      Long,               Lat,       h
    [x1,y1,z1]= sph2cart(deg2rad(Point1(:,2)), deg2rad(Point1(:,1)), Re);
    [x2,y2,z2]= sph2cart(deg2rad(Point2(2)), deg2rad(Point2(1)), Re);
    CartDistance = sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
    CartDistance = CartDistance./1000; % [km]   
end 