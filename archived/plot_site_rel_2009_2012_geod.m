function plot_site_rel_2009_2012_geod(site,median_smoothing_window,gaussian_smoothing_window,iterations,sigma_mult,athome);

%% Original code written by Matt King, Newcastle University. Modifications by Ian Bartholomew and Andrew Sole, University of Edinburgh.

% Modified by Andrew Tedstone October 2012. 
% - Flexible time buffer, duration worked out from data rather than hard coded.
% - Deletes Day 116 2012
% - Options (commented out) for dealing with special cases of lev7 and lev8.


% useage:
% plot_site_rel_2009_2012_geod('lev4',7200,7200,2,2,0);
%site='lev4';
%median_smoothing_window=7200;
%gaussian_smoothing_window=7200;
%iterations=2;
%sigma_mult=2;
%athome=0;

% site = gps site name (e.g. 'lev6');
% median_smoothing_interval = smoothing window length for median filter in seconds
% gaussian_smoothing_interval = smoothing window length for gaussian filter in seconds
% iterations = number of times to run median filter
% sigma_mult = data points beyond sigma_mult*sigma will be removed during smoothing
% athome = whether or not to use signal processing toolbox functions (i.e. dont use at work (athome=0) to avoid licensing issues)

% load existing velocity data
disp('loading 2009 velocity')
load ../gps_2009_2010_velocities/2009_velocity
% choose relevant site data
switch lower(site)
   case {'lev0'}
      vel_2009 = v1(:,1:2);
	  load ../gps_2009_2010_velocities/lev0/lev0_velocity_2010;
	  vel_2010 = v_24h;
   case {'lev1'}
      vel_2009 = v2(:,1:2);
	  load ../gps_2009_2010_velocities/lev1/lev1_velocity_2010;
	  vel_2010 = v_24h; 
   case {'lev2'}
      vel_2009 = v3(:,1:2);
	  load ../gps_2009_2010_velocities/lev2/lev2_velocity_2010;
	  vel_2010 = v_24h;	  
   case {'lev3'}
      vel_2009 = v4(:,1:2);
	  load ../gps_2009_2010_velocities/lev3/lev3_velocity_2010;
	  vel_2010 = v_24h;	  
   case {'lev4'}
      vel_2009 = v5(:,1:2);
	  load ../gps_2009_2010_velocities/lev4/lev4_velocity_2010;
	  vel_2010 = v_24h;	  
   case {'lev5'}
      vel_2009 = v6(:,1:2);
	  load ../gps_2009_2010_velocities/lev5/lev5_velocity_2010;
	  vel_2010 = v_24h;	  
   case {'lev6'}
      vel_2009 = v7(:,1:2); 
	  load ../gps_2009_2010_velocities/lev6/lev6_velocity_2010;
	  vel_2010 = v_24h;     
end
% add on 365 days to doy for 2010
%vel_2010(:,1) = vel_2010(:,1)+365;
% 
% clear unnecessary variables
clear v1* v2* v3* v4* v5* v6* v7* v_6h v_24h i_loop

%% set variables
file=([site '_2012_geod.dat']);

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
smap(:,18) = abs(smap(:,18));
smap(:,19) = abs(smap(:,19));

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

%% Compute 2d velocity magnitude and covariance
% North velocity (assuming linear)
[b,Cb,vf]=lreg(smap(1:end,14)./365.25,smap(1:end,18));
vel_n=b(2);

%% Now rotate into along- and across-track displacements
% find direction
fprintf('Rotating into along- and cross-track pairs\n');
dir=atan2(vel_n,vel_e);
%rotate --> xy(:,1)=along track, xy(:,2)=across track.
R1=[cos(-dir) -sin(-dir); sin(-dir) cos(-dir)];
xy=(R1*[smap(:,19) smap(:,18)]')';

%% Data smoothing
% find data interval. use mode to account for data gaps due to removal of bad values
for j=1:100
	interval(j)=abs(smap(j+1,3)-smap(j,3));
end
interval = mode(interval);
fprintf(['Interval = ',num2str(interval),' seconds\n']);
% converts doy column into date vector
fprintf('Managing timeseries data\n');
tvec=datevec(smap(:,14));
tvec(:,6)=round(tvec(:,6));
% returns date vector to doy at matlab resolution 
tnum=round2(datenum(tvec),1e-8);
% make ideal time series intervals
ti=(datenum(tvec(1,:)):(datenum(0,1,0,0,0,interval)):datenum(tvec(end,:)));
ti_vec=datevec(ti);
ti_vec(:,6)=round(ti_vec(:,6));
ti=round2(datenum(ti_vec),1e-8);

% Interpolate missing values in time series
% get median y, z and time values for data points with identical x values as interp1 needs distinct x values.
disp('Consolidating data based on x values');
[xcon,ycon,ind] = consolidator(xy(:,1),xy(:,2),'median');
[xcon,zcon,ind] = consolidator(xy(:,1),smap(:,20),'median');
[xcon,tcon1,ind] = consolidator(xy(:,1),tnum,'median');
% get median x, y and z values for data points with identical time values as interp1 needs distinct x values.
disp('Consolidating data based on time values');
[tcon,xcon,ind] = consolidator(tcon1,xcon,'median');
[tcon,ycon,ind] = consolidator(tcon1,ycon,'median');
[tcon,zcon,ind] = consolidator(tcon1,zcon,'median');

% need to do initial removal of bad data before interpolating, otherwise some interpolated points are also bad, meaning that the subsequent median filtering is less effective.

xcon = mediannan(xcon,15);
ycon = mediannan(ycon,15);
zcon = mediannan(zcon,15);

%nano = find(isnan(xcon) | isnan(ycon) | isnan(zcon));
%if isempty(nano) ~= 1
%	zcon2filt = zcon(~nano);
%	zcon_filt = RankOrderFilter(zcon2filt,240,50);
%	zcon(~nano) = zcon_filt;
%else
%	zcon = RankOrderFilter(zcon,240,50);
%end



% interpolate consolidated values
xi=interp1(tcon,xcon,ti,'linear'); % xi is now along track displacement
yi=interp1(tcon,ycon,ti,'linear'); % yi is now across track displacement
zi=interp1(tcon,zcon,ti,'linear');
% get rid of any nans...not sure why there are any, but sometimes there are?
nano = find(isnan(xi));
xi(nano) = [];
yi(nano) = [];
zi(nano) = [];
ti(nano) = [];
nano = find(isnan(yi));
xi(nano) = [];
yi(nano) = [];
zi(nano) = [];
ti(nano) = [];
nano = find(isnan(zi));
xi(nano) = [];
yi(nano) = [];
zi(nano) = [];
ti(nano) = [];
nano = find(isnan(ti));
xi(nano) = [];
yi(nano) = [];
zi(nano) = [];
ti(nano) = [];
%
[b,Cb,vf]=lreg(ti,xi);
vel_x=b(2);
xi_detrended=xi-(ti*b(2)+b(1));

%% Remove outliers
fprintf('Removing outliers\n');
for r=1:iterations 
    if r==1
        fprintf(['Iteration ',num2str(r),' of ' num2str(iterations)]);
        xs=RankOrderFilter(xi_detrended(1:(990/interval):end),(990/interval),50);  % 50 = 50th percentile = median
        xs=interp1(ti(1:(990/interval):end),xs,ti,'linear');
        fprintf('. ');
        ys=RankOrderFilter(yi(1:(990/interval):end),(990/interval),50);
        ys=interp1(ti(1:(990/interval):end),ys,ti,'linear');
        fprintf('. ');
        zs=RankOrderFilter(zi(1:(990/interval):end),(990/interval),50);
        zs=interp1(ti(1:(990/interval):end),zs,ti,'linear');
        fprintf('. ');
        xres=xi_detrended-xs;
        yres=yi-ys;
        zres=zi-zs;
        i=find((abs(xres))>=0.08|(abs(yres))>=0.04|(abs(zres))>=0.15); % what are these values based on?
    else
		fprintf(['Iteration ',num2str(r),' of ' num2str(iterations)]);
		xs=RankOrderFilter(xi_detrended,median_smoothing_window/interval,50);
		fprintf('. ');
		ys=RankOrderFilter(yi,median_smoothing_window/interval,50);
		fprintf('. ');
		zs=RankOrderFilter(zi,median_smoothing_window/interval,50);
		fprintf('. ');
		xres=xi_detrended-xs;
		yres=yi-ys;
		zres=zi-zs;
		sig_x=movingstd(xres,median_smoothing_window/interval);
		fprintf('. ');
		sig_y=movingstd(yres,median_smoothing_window/interval);
		fprintf('. ');
		sig_z=movingstd(zres,median_smoothing_window/interval);
		fprintf('. ');
		i=find((abs(xres))>=(sig_x*sigma_mult)|(abs(yres))>=(sig_y*sigma_mult)|(abs(zres))>=(sig_z*sigma_mult));
    end
    tj=ti;
    tj(i)=[];
    xi_detrended(i)=[];
    xi_detrended=interp1(tj,xi_detrended,ti,'linear');
    yi(i)=[];
    yi=interp1(tj,yi,ti,'linear');
    zi(i)=[];
    zi=interp1(tj,zi,ti,'linear');
    fprintf('done\n');
end
xi=xi_detrended+(ti*b(2)+b(1));

%% now can do the smoothing using gaussian lowpass filter
if athome == 1
	fprintf('Smoothing data\n');
	xf=filtfilt((gaussfiltcoef((1/interval),(1/(gaussian_smoothing_window)))),1,xi);
	yf=filtfilt((gaussfiltcoef((1/interval),(1/(gaussian_smoothing_window)))),1,yi);
	zf=filtfilt((gaussfiltcoef((1/interval),(1/(gaussian_smoothing_window*4)))),1,zi);
	xf_detrended=xf-(ti*b(2)+b(1));
else
	xf=xi;
	yf=yi;
	zf=zi;
	xf_detrended=xf-(ti*b(2)+b(1));
end

% remove interpolated values
fprintf('Removing interpolated values\n');
i=find(ismember(ti,tnum));
xyz=[xf(i), yf(i), zf(i)];
xyzt=[xf(i), yf(i), zf(i), ti(i)];

% plot various stages of filtering
figure
plot(xy(:,1),xy(:,2),'.k');
hold on
plot(xcon,ycon,'.');
plot(xi,yi,'.r');
plot(xyz(:,1),xyz(:,2),'.c');
legend('raw','initial filtering + consolidated','interpolated','filtered');
xlabel('X (m)');
ylabel('Y (m)');
title('X vs. Y');
axis equal

% plot various stages of filtering of up/z
figure
plot(xy(:,1),smap(:,20),'.k');
hold on
plot(xcon,zcon,'.');
plot(xi,zi,'.r');
plot(xyz(:,1),xyz(:,3),'.c');
legend('raw','initial filtering + consolidated','interpolated','filtered');
xlabel('X (m)');
ylabel('Z (m)');
title('X vs. Up');
axis equal

%% Take linear fit from z-data
[b,Cb,vf]=lreg(xyz(:,1),xyz(:,3));
z_fit(:,1)=ti;
z_fit(:,2)=zf-xf.*b(2);

%% Calculate mean daily velocities by differencing
fprintf('Computing daily velocities\n');
 i=find(ti==ceil(ti(1)));
v_24h(:,1)=ti(i:(86400/interval):end-(86400/interval));
v_24h(:,2)=(xf(i+(86400/interval):(86400/interval):end)-xf(i:(86400/interval):end-(86400/interval)))*365;

%% buffer edges so that dataset runs through whole-year lengths of d/s
% AJT addition 23 oct 2012
start_year = int16(smap(1,1))
start_year = 2009;
end_year = int16(smap(end,1))
buffer_l = 0;
for yr=start_year:end_year
  if eomday(yr,2) == 29
    buffer_l = buffer_l + 366;
  else
    buffer_l = buffer_l + 365;
  end
end
buffer_l
buffer(:,1) = (1:buffer_l);
buffer(:,2) = buffer(:,1)*NaN;
% insert new velocity data
buffer(v_24h(:,1),2) = v_24h(:,2);
% create new buffered dataset
v_24h = buffer;

%% create an alternative 2009-2012 dataset using existing velocity data where available
% get days where velocity data already exists
% insert = ([vel_2009(:,1);vel_2010(:,1)]);
% % create new variable
% v_24h_combo = v_24h;
% % insert existiing data into correct days of new variable
% v_24h_combo(insert,2) = ([vel_2009(:,2);vel_2010(:,2)]);
% % find and replace nans with data from new dataset which does not have nans
% nano = find(isnan(v_24h_combo(:,2)));
% v_24h_combo(nano,2) = v_24h(nano,2);

% %% make a combination dataset with existing 2009 data and new 2010-2011 data
% v_24h_opt = v_24h;
% idx2009 = find(v_24h_opt(:,1)<=365);
% v_24h_opt(idx2009,2) = v_24h_combo(idx2009,2);

% get rid of strange value on day 517 in new 2010 data
 switch lower(site)
    case {'lev1'}
 		v_24h_opt(517,2) = v_24h_combo(517,2);
 end

% get rid of data where lev0 redrilled without measuring position difference
switch lower(site)
   case {'lev0'}
		v_24h(930:974,:) = NaN;
		v_24h_opt(930:974,:) = NaN;
		v_24h_combo(930:974,:) = NaN;
	case {'lev1'}
	  %v_24h(930:975,:) = NaN;
  	%v_24h_opt(930:975,:) = NaN;
	  %v_24h_combo(930:975,:) = NaN;		
end

% get rid of 'incorrect' lev6 velocity on aug 2012 redrill day
if site == 'lev6'
  v_24h(1338,2) = NaN;
  v_24h_opt(1338,2) = NaN;
  v_24h_combo(1338,2) = NaN;
end

%%  Calcualate short term velocities by differencing
fprintf('Computing short term velocities\n');
fprintf('Differencing\n');
v_6h=zeros(length(xf)-(21600/interval),3);
for i=1:(length(xf)-(21600/interval))
v_6h(i,1)=ti(i+(21600/(2*interval)));
v_6h(i,2)=(xf(i+21600/interval)-xf(i))*1460;
end
v_4h=zeros(length(xf)-(14400/interval),3);
for i=1:(length(xf)-(14400/interval))
v_4h(i,1)=ti(i+(14400/(2*interval)));
v_4h(i,2)=(xf(i+14400/interval)-xf(i))*2190;
end

%% save velocity data
savefile = ([site 'vel_2009_spr2013_UTC-2_20130607.mat']);
if (site == 'lev7') | (site == 'lev8')
  save(savefile, 'v_24h','v_6h')
else
  save(savefile, 'v_24h','v_6h','v_24h_opt','v_24h_combo')
end
%% save xyzt data
savefile = ([site 'xyzt_2009_spr2013_UTC-2_20130607.mat']);
save(savefile, 'xyzt')

fprintf('Finished\n');

%% plotting
% set ticks and tick labels
allxtick = [1 32 60 91 121 152 182 213 244 274 305 335 366 397 425 456 486 517 547 578 609 639 670 700 731 762 790 821 851 882 912 943 974 1004 1035 1065 1096 1127 1156 1187];
allxticklabels = {'1/2009','2/2009','3/2009','4/2009','5/2009','6/2009','7/2009','8/2009' ...
,'9/2009','10/2009','11/2009','12/2009','1/2010','2/2010','3/2010','4/2010','5/2010','6/2010','7/2010' ...
,'8/2010','9/2010','10/2010','11/2010','12/2010','1/2011','2/2011','3/2011','4/2011','5/2011','6/2011','7/2011' ...
,'8/2011','9/2011','10/2011','11/2011','12/2011','1/2012','2/2012','3/2012','4/2012','5/2012'};
ticks2use = 4;
% velocity and height
figure
% get rid on nans
zf_plot = xyzt(:,3);
ti_plot = xyzt(:,4);
nano = isnan(xyzt(:,3));
zf_plot(nano,:) = [];
ti_plot(nano,:) = [];
grid on
[AX1, H11, H12]=plotyy(v_24h(:,1),v_24h(:,2),ti_plot,zf_plot,'stairs','plot');
set(AX1(1),'xlim',[1 1220],'xtick',[1 32 60 91 121 152 182 213 244 274 305 335 366 397 425 456 486 517 547 ...
578 609 639 670 700 731 762 790 821 851 882 912 943 974 1004 1035 1065 1096 1127 1156 1187], ...
'xticklabel',{'1/2009','2/2009','3/2009','4/2009','5/2009','6/2009','7/2009','8/2009' ...
,'9/2009','10/2009','11/2009','12/2009','1/2010','2/2010','3/2010','4/2010','5/2010','6/2010','7/2010' ...
,'8/2010','9/2010','10/2010','11/2010','12/2010','1/2011','2/2011','3/2011','4/2011','5/2011','6/2011','7/2011' ...
,'8/2011','9/2011','10/2011','11/2011','12/2011','1/2012','2/2012','3/2012','4/2012','5/2012'},'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'ycolor','b')
set(H11,'Color','b','linewidth',1)
ylabel(AX1(1),{'velocity (my^-^1)'},'fontsize',10)
title(site,'fontsize',10)
set(AX1(2),'xlim',[1 1220],'xtick',[1 32 60 91 121 152 182 213 244 274 305 335 366 397 425 456 486 517 547 ...
578 609 639 670 700 731 762 790 821 851 882 912 943 974 1004 1035 1065 1096 1127 1156 1187], ...
'xticklabel',{'1/2009','2/2009','3/2009','4/2009','5/2009','6/2009','7/2009','8/2009' ...
,'9/2009','10/2009','11/2009','12/2009','1/2010','2/2010','3/2010','4/2010','5/2010','6/2010','7/2010' ...
,'8/2010','9/2010','10/2010','11/2010','12/2010','1/2011','2/2011','3/2011','4/2011','5/2011','6/2011','7/2011' ...
,'8/2011','9/2011','10/2011','11/2011','12/2011','1/2012','2/2012','3/2012','4/2012','5/2012'},'fontsize',10,'ycolor','k')
ylabel(AX1(2),{'surface';' height (m)'},'fontsize',10)
set(H12,'Color','k','LineStyle','none','Marker','o','MarkerSize',2)
linkaxes([AX1(1) AX1(2)],'x');
% Print figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
% A4 = 8.3 x 11.7 inches
set(gcf, 'PaperPosition', [0.25 0.25 11.2 7.8]);
% save
print -depsc2 test.eps
fixPSlinestyle('test.eps',[site '_vel_elev_2009-2013.eps']);

% different velocity options
figure;
subplot(3,1,2);
stairs(v_24h(:,1),v_24h(:,2),'linewidth',2,'color',[.7 .7 .7])
hold on                                                                    
stairs(v_24h_combo(:,1),v_24h_combo(:,2),'k'); 
stairs(v_24h_opt(:,1),v_24h_opt(:,2),'k');
title(site);
ylabel(gca,{'velocity (my^-^1)'},'fontsize',10)
% set x and y ticks
switch lower(site)
   case {'lev0'}
        set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'xlim',[1 1220]);
   case {'lev1'}
        set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'xlim',[1 1220]);
   case {'lev2'}
        set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'xlim',[1 1220]);
   case {'lev3'}
        set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'xlim',[1 1220]);
   case {'lev4'}
        set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'xlim',[1 1220]);
   case {'lev5'}
        set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'xlim',[1 1220]);
   case {'lev6'}
       set(gca,'xtick',allxtick([1:ticks2use:end]),'xticklabel',allxticklabels([1:ticks2use:end]),'fontsize',10,'ylim',[0 100],'ytick',[0 20 40 60 80 100],'xlim',[1 1220]);
end
set(gca,'xlim',[121 1220]);
% add background velocity
switch lower(site)
   case {'lev0'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[48.8 48.8],'color',[0.7 0.7 0.7]);
   case {'lev1'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[102.5 102.5],'color',[0.7 0.7 0.7]);	
   case {'lev2'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[94.8 94.8],'color',[0.7 0.7 0.7]);
   case {'lev3'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[114.2 114.2],'color',[0.7 0.7 0.7]);
   case {'lev4'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[108.5 108.5],'color',[0.7 0.7 0.7]);
   case {'lev5'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[88.5 88.5],'color',[0.7 0.7 0.7]);
   case {'lev6'}
		plot([min(get(gca,'xlim')) max(get(gca,'xlim'))],[61.9 61.9],'color',[0.7 0.7 0.7]);	
end
%legend('new','old with new gap-fill','old 2009 new 2010');
% Print figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
% A4 = 8.3 x 11.7 inches
set(gcf, 'PaperPosition', [0.25 0.25 11.2 7.8]);
% save
print -depsc2 test.eps
fixPSlinestyle('test.eps',[site '_vel_options_2009-2013.eps']);