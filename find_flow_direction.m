site = 'lev7';

approach = 'existing'; % or test

%% set variables
file=([site '_2011_2013_geod.dat']);

%% load the data
% trailing \n removed as python doesn't output this linebreak.
fprintf(['Loading data file: ',file,'\n']);
[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20]=textread(file,'%d %d %f %f %f %f %f %f %f %f %d %f %f %f %d %d %d %f %f %f'); 
smap=[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20];
clear s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s20

%% fix years if multiple years in dataset (only if not already done somewhere else)
%smap(:,14)=smap(:,14)+(smap(:,1)-min(smap(:,1)))*365;

%% For Lev7 2011-2012 ONLY - pad up fractDOY to make it the same as all other datasets
if site == 'lev7'
  smap(:,14) = smap(:,14) + (365 * 2);
end

%% For Lev8 2012-2013 ONLY - pad up fractDOY to make it the same as all other datasets
if site == 'lev8'
  smap(:,14) = smap(:,14) + (365 * 3);
end


%% Remove 2012 redrill days - shouldn't be necessary but they cause short-duration spikes
% This only works for sites which continue recording beyond this day - so e.g. NOT lev6 autumn 2012.
disp('Removing day 1211 (116 2012)')
get_rid = find(smap(:,14) >= 1211 & smap(:,14) < 1212);
smap(get_rid,:) = [];

if site ~= 'lev6'
  disp('Removing day 1338 (243 2012)')
  get_rid = find(smap(:,14) >= 1338 & smap(:,14) < 1339);
  smap(get_rid,:) = [];
end


%% Remove site-specific dates
if site == 'lev0'
  disp('Removing various from lev0 winter-spring 2013')
  get_rid = find(smap(:,14) >= 1486 & smap(:,14) < 1487.9);
  smap(get_rid,:) = [];
  get_rid = find(smap(:,14) >= 1493.7 & smap(:,14) < 1498);
  smap(get_rid,:) = [];
  get_rid = find(smap(:,14) >= 1504 & smap(:,14) < 1505.7);
  smap(get_rid,:) = [];
  get_rid = find(smap(:,14) >= 1512.4 & smap(:,14) < 1513.71);
  smap(get_rid,:) = [];
  get_rid = find(smap(:,14) >= 1520.4 & smap(:,14) < 1522.34);
  smap(get_rid,:) = [];
  get_rid = find(smap(:,14) >= 1530.5 & smap(:,14) < 1531.6);
  smap(get_rid,:) = [];
  get_rid = find(smap(:,14) >= 1533.7 & smap(:,14) < 1535.98);
  smap(get_rid,:) = [];
end

if site == 'lev1'
  disp('Removing day 1375-1377 lev1')
  get_rid = find(smap(:,14) >= 1375 & smap(:,14) < 1377);
  smap(get_rid,:) = [];
end

if site == 'lev2'
  disp('Removing day 1417-1419 lev2')
  get_rid = find(smap(:,14) >= 1416 & smap(:,14) < 1420);
  smap(get_rid,:) = [];
end

if site == 'lev5'
  disp('Removing day 1385-1387 lev5')
  get_rid = find(smap(:,14) >= 1385 & smap(:,14) < 1388);
  smap(get_rid,:) = [];
end

if site == 'lev7'
  disp('Removing day 239 2011 Lev7')
  get_rid = find(smap(:,14) >= (239+365+365) & smap(:,14) < (240+365+365));
  smap(get_rid,:) = [];
end

if site == 'lev8'
  disp('Removing day 1347-1349 lev8')
  get_rid = find(smap(:,14) >= 1347 & smap(:,14) < 1350);
  smap(get_rid,:) = [];
  disp('Removing day 1393-1398 lev8')
  get_rid = find(smap(:,14) >= 1393 & smap(:,14) < 1399);
  smap(get_rid,:) = [];
  disp('Removing day 1407 lev8')
  get_rid = find(smap(:,14) >= 1407 & smap(:,14) < 1408);
  smap(get_rid,:) = [];
  disp('Removing day 1415 lev8')
  get_rid = find(smap(:,14) >= 1415 & smap(:,14) < 1416);
  smap(get_rid,:) = [];
  disp('Removing day 1422-1425 lev8')
  get_rid = find(smap(:,14) >= 1422 & smap(:,14) < 1425);
  smap(get_rid,:) = [];
end



% sort based on day of year (only if not already sorted);
%smap = sortrows(smap,14);

% make displacements absolute
smap_nonabs = smap;
%smap(:,18) = abs(smap(:,18));
%smap(:,19) = abs(smap(:,19));


% remove bad data where RMS is large (>60mm) or height std dev is high
% (>10cm), or RMS is NaN
i=find((smap(:,10)>90)|(smap(:,9)>90)|(isnan(smap(:,10))));
fprintf('Editing out %d bad epochs\n',size(i,1));
smap(i,:)=[];

% remove epochs where filter has interpolated values (generally due to wrong sampling interval set)
i=find(smap(:,11)==0);
fprintf('Editing out %d epochs with no double differences (normally due to wrong sampling interval settings)\n',size(i,1));
smap(i,:)=[];

%% for Greenland, set time to local (UTC-2)
fprintf('Setting local time in GPS data\n');
smap(:,14)=smap(:,14)-(2/24);

% scale down to avoid inversion of large matrices
% if (size(smap(:,14))>20000)
%     step = 1000;
% elseif (size(smap(:,14))>2000)
%     step = 100;
% else
%     step = 1;
% end

%% Compute 2d velocity magnitude and covariance
% East velocity (assuming linear)

[b,Cb,vf]=lreg(smap(1:end,14)./365.25,smap(1:end,19));
vel_e=b(2);
intercept_e = b(1);

%% Compute 2d velocity magnitude and covariance
% North velocity (assuming linear)

[b,Cb,vf]=lreg(smap(1:end,14)./365.25,smap(1:end,18));
vel_n=b(2);
intercept_n = b(1);

%% Now rotate into along- and across-track displacements
% find direction
fprintf('Rotating into along- and cross-track pairs\n');
dir=atan2(vel_n,vel_e);
%rotate --> xy(:,1)=along track, xy(:,2)=across track.
%R1=[cos(dir) -sin(dir); sin(dir) cos(dir)];
R1=[cos(-dir) -sin(-dir); sin(-dir) cos(-dir)];
xy=(R1*[smap(:,19) smap(:,18)]')';


% ALTERNATIVE: find direction using linear regression 
% between east and north coordinates
[b,Cb,vf] = lreg(smap(1:end,19),smap(1:end,18));
slope = b(2);
intercept = b(1);

figure()
% Plot original east north
plot(smap(:,19),smap(:,18),'xb')
hold on
% plot direction from original method
x = [min(smap(:,19)):1:max(smap(:,19))];
y = (vel_n / vel_e) * x;
plot(x,y,'g')
% Plot direction using E-N regression
y1 = slope * x;
plot(x,y1,'r')

% Plot rotated data from original method
plot(xy(:,1),xy(:,2),'xk');
% Plot a zero y line
plot(x,0,'g');



