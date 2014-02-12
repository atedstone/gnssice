% gps_positions_heli

% pete's estimates of leverett gps sites
% site 2 (Lev 1)
lat_p(:,1) = [67 5.576];
lon_p(:,1) = [50 1.766];
% site 3
lat_p(:,2) = [67 6.160];
lon_p(:,2) = [49 48.474];
% site 4
lat_p(:,3) = [67 6.917];
lon_p(:,3) = [49 28.039];
% site 5
lat_p(:,4) = [67 7.583];
lon_p(:,4) = [49 0.850];
% site 6
lat_p(:,5) = [67 9.183];
lon_p(:,5) = [48 22.433];
% site 7
lat_p(:,6) = [67 9.755];
lon_p(:,6) = [47 33.067];
% site 8
lat_p(:,7) = [67 8.421];
lon_p(:,7) = [49 48.408];
% site 9 (lev8)
lat_p(:,8) = [67 7.345];
lon_p(:,8) = [49 48.462];

% from processed rinex files for end of 2010 season
% site 2 (lev 1)
lat(:,1) = [67 5.5675];
lon(:,1) = [50 1.8343];
% site 3
lat(:,2) = [67 6.2612];
lon(:,2) = [49 48.6680];
% site 4
lat(:,3) = [67 6.9234];
lon(:,3) = [49 24.2788];
% site 5
lat(:,4) = [67 7.5817];
lon(:,4) = [49 0.8419];
% site 6
lat(:,5) = [67 9.1903];
lon(:,5) = [48 22.4245];
% site 7
lat(:,6) = [67 9.7563];
lon(:,6) = [47 33.0624];

% from processed rinex files for end of 2011 season
% site	year	day		lat				lon				elevation
% lev0	2011	231		67.06957542  	309.86979057   	446.518
% lev1	2011	207		67.09268787  	309.96683091   	607.909
% lev2	2011	238		67.10418598  	310.18648243   	791.349
% lev3	2011	243		67.11515407  	310.59250503  	1054.856
% lev4	2011	244		67.12617104  	310.98315585  	1217.939
% lev5	2011	244		67.15329788  	311.62407887  	1475.145
% lev6	2011	244		67.16263899  	312.44745134  	1714.651
% lev7	2011	240		67.14634979  	310.18749251   	831.805
% lev8	-		-		67.1224			310.1923		-

% lev0
lat_2011(:,1) = 67.06957542;
lon_2011(:,1) = 360-309.86979057;
% lev1
lat_2011(:,2) = 67.09268787;
lon_2011(:,2) = 360-309.96683091;
% lev2
lat_2011(:,3) = 67.10418598;
lon_2011(:,3) = 360-310.18648243;
% lev3
lat_2011(:,4) = 67.11515407;
lon_2011(:,4) = 360-310.59250503;
% lev4
lat_2011(:,5) = 67.12617104;
lon_2011(:,5) = 360-310.98315585;
% lev5
lat_2011(:,6) = 67.15329788;
lon_2011(:,6) = 360-311.62407887;
% lev6
lat_2011(:,7) = 67.16263899;
lon_2011(:,7) = 360-312.44745134;
% lev7
lat_2011(:,8) = 67.14634979;
lon_2011(:,8) = 360-310.18749251;
% lev8
lat_2011(:,9) = 67.1224;
lon_2011(:,9) = 360-310.1923;


% for spring 2013, from end of rinex files for 2012 season
% lev0
lat_2012(:,1) = 67.06925201;
lon_2012(:,1) = 360-309.868888615;
% lev1
lat_2012(:,2) = 67.092595972;
lon_2012(:,2) = 360-309.964172452;
% lev2
lat_2012(:,3) = 67.104085211;
lon_2012(:,3) = 360-310.184111881;
% lev3
lat_2012(:,4) = 67.114935149;
lon_2012(:,4) = 360-310.589812615;
% lev4
lat_2012(:,5) = 67.126005092;
lon_2012(:,5) = 360-310.98061344;
% lev5
lat_2012(:,6) = 67.153407634;
lon_2012(:,6) = 360-311.62197273;
% lev6
lat_2012(:,7) = 67.162663555;
lon_2012(:,7) = 360-312.44598445;
% lev7
lat_2012(:,8) = 67.14620005;
lon_2012(:,8) = 360-310.185888888;
% lev8
lat_2012(:,9) = 67.128863764;
lon_2012(:,9) = 360-310.184392065;



%% convert ddm to dd
lat_p_dd = dm2degrees(lat_p');
lon_p_dd = dm2degrees(lon_p');

lat_dd = dm2degrees(lat');
lon_dd = dm2degrees(lon');

lat_2011_dm = degrees2dm(lat_2011(:));
lon_2011_dm = degrees2dm(lon_2011(:));

lat_2012_dm = degrees2dm(lat_2012(:));
lon_2012_dm = degrees2dm(lon_2012(:));

%% plot

figure;
plot(-lon_p_dd,lat_p_dd,'+k');
hold on
plot(-lon_dd,lat_dd,'or');
plot(-lon_2011,lat_2011,'sb');
plot(-lon_2012,lat_2012,'^b');
%plot(-lon_2011(end),lat_2011(end),'xb');
% plot lev8
%lev8_lat = dm2degrees([67 07.738])
%lev8_lon = dm2degrees([49 48.879])
%plot(-lev8_lon,lev8_lat,'^');
title('Geodetic');
legend('PWN estimates','rinex end 2010','rinex end 2011','rinex end 2012');

%% convert to utm
[x_p,y_p,utmzone_p] = deg2utm(lat_p_dd,-lon_p_dd);
[x,y,utmzone] = deg2utm(lat_dd,-lon_dd);
[x_2011,y_2011,utmzone] = deg2utm(lat_2011,-lon_2011);
[x_2012,y_2012,utmzone] = deg2utm(lat_2012,-lon_2012);

figure;
plot(x_p,y_p,'+k');
hold on
plot(x,y,'or');
plot(x_2011,y_2011,'sb');
plot(x_2011(end),y_2011(end),'xb');
plot(x_2012(end),y_2012(end),'^b');
title('UTM');
legend('PWN estimates','rinex end 2010','rinex end 2011','rinex end 2012');

%% calculate distance between points
for i = 1:6
    dist_diff(i) = sqrt((x_p(i)-x(i)).^2+(y_p(i)-y(i)).^2);
end
dist_diff(:)/1000
