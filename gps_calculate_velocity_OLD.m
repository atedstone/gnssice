%% Original code written by Matt King, Newcastle University. 
% Modifications by Ian Bartholomew and Andrew Sole, University of Edinburgh.
% Modified by Andrew Tedstone, University of Edinburgh, September 2012.

% gps_calculate_velocity.m (was plot_site_rel_2009_2011_geod.m)

clear

site = 'lev2';
year = '2009-2012'
%site = input('Site: ','s');
%year = input('Year: ','s');

%% set variables
file=([site '_' year '_GEOD_m.dat']);

%% load the data
fprintf(['Loading data file: ',file,'\n']);
[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20]=textread(file,'%d %d %f %f %f %f %f %f %f %f %d %f %f %f %d %d %d %f %f %f'); 
smap=[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20];
clear s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17 s18 s19 s20

%% fix years if multiple years in dataset (only if not already done somewhere else)
%smap(:,14)=smap(:,14)+(smap(:,1)-min(smap(:,1)))*365;

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
j=find(smap(:,11)==0);
fprintf('Editing out %d epochs with no double differences (normally due to wrong sampling interval settings)\n',size(j,1));
smap(j,:)=[];

%k = find(smap(:,19) < 0);
%smap(k,:) = [];

%% for Greenland, set time to local (UTC-2) D
%fprintf('Setting local time in GPS data\n');
%smap(:,14)=smap(:,14)-(2/24);

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

% get rid of bad lev0 data
%if site == 'lev0'
%	get_rid = find(xy(:,2)>10);
%	xy(get_rid,:) = [];
%	smap(get_rid,:) = [];
%end

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

%% Interpolate missing values in time series
% get median y, z and time values for data points with identical x values as interp1 needs distinct x values.
[xcon,ycon,ind] = consolidator(xy(:,1),xy(:,2),'median');
[xcon,zcon,ind] = consolidator(xy(:,1),smap(:,20),'median');
[xcon,tcon1,ind] = consolidator(xy(:,1),tnum,'median');
disp('Consolidating data based on x values');
% get median x, y and z values for data points with identical time values as interp1 needs distinct x values.
[tcon,xcon,ind] = consolidator(tcon1,xcon,'median');
[tcon,ycon,ind] = consolidator(tcon1,ycon,'median');
[tcon,zcon,ind] = consolidator(tcon1,zcon,'median');
disp('Consolidating data based on time values');

% need to do initial removal of bad data before interpolating, otherwise some interpolated points are also bad, meaning that the subsequent median filtering is less effective.
xcon = mediannan(xcon,15);
ycon = mediannan(ycon,15);
zcon = mediannan(zcon,15);

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
num=2; % increase to smooth more times. Ian used num=2
for r=1:num 
    if r==1
        fprintf(['Iteration ',num2str(r),' of ' num2str(num)]);
        xs=RankOrderFilter(xi_detrended(1:(990/interval):end),(990/interval),50);
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
        i=find((abs(xres))>=0.08|(abs(yres))>=0.04|(abs(zres))>=0.15);
    else
		fprintf(['Iteration ',num2str(r),' of ' num2str(num)]);
		xs=RankOrderFilter(xi_detrended,7200/interval,50);
		fprintf('. ');
		ys=RankOrderFilter(yi,7200/interval,50);
		fprintf('. ');
		zs=RankOrderFilter(zi,7200/interval,50);
		fprintf('. ');
		xres=xi_detrended-xs;
		yres=yi-ys;
		zres=zi-zs;
		sig_x=movingstd(xres,7200/interval);
		fprintf('. ');
		sig_y=movingstd(yres,7200/interval);
		fprintf('. ');
		sig_z=movingstd(zres,7200/interval);
		fprintf('. ');
		mult=2; % reduce to remove more outliers. Ian used mult=2
		i=find((abs(xres))>=(sig_x*mult)|(abs(yres))>=(sig_y*mult)|(abs(zres))>=(sig_z*mult));
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
% N.b. in Andrew's script this doesn;t happen if being run in dept (licensing)
% fprintf('Smoothing data\n');
% xf=filtfilt((gaussfiltcoef((1/interval),(1/(7200*1)))),1,xi);
% yf=filtfilt((gaussfiltcoef((1/interval),(1/(7200*1)))),1,yi);
% zf=filtfilt((gaussfiltcoef((1/interval),(1/(7200*4)))),1,zi);
% xf_detrended=xf-(ti*b(2)+b(1));

	xf=xi;
	yf=yi;
	zf=zi;
	xf_detrended=xf-(ti*b(2)+b(1));

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
legend('raw','consolidated','interpolated','filtered');
xlabel('X (m)');
ylabel('Y (m)');
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
start_year = int16(smap(1,1))
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

%%  Calcualate short term velocities by differencing
fprintf('Computing short term (6h) velocities\n');
fprintf('Differencing (using sliding window)\n');

v_6h=zeros(length(xf)-(21600/interval),3);
for i=1:(length(xf)-(21600/interval))
  v_6h(i,1)=ti(i+(21600/(2*interval)));
  v_6h(i,2)=(xf(i+21600/interval)-xf(i))*1460;
end
% v_4h=zeros(length(xf)-(14400/interval),3);
% for i=1:(length(xf)-(14400/interval))
% v_4h(i,1)=ti(i+(14400/(2*interval)));
% v_4h(i,2)=(xf(i+14400/interval)-xf(i))*2190;
% end


%% save velocity data
savefile = ([site '_vel_' year '-myscript-r2.mat']);
save(savefile, 'v_24h','v_6h')
%% save xyzt data
savefile = ([site '_xyzt_' year '-myscript-r2.mat']);
save(savefile, 'xyzt')

fprintf('Finished.\n');

%% plotting
% set ticks and tick labels
% allxtick = [1 32 60 91 121 152 182 213 244 274 305 335];
% allxticklabels = {'1/2011','2/2011','3/2011','4/2011','5/2011','6/2011','7/2011' ...
% ,'8/2011','9/2011','10/2011','11/2011','12/2011'};
% ticks2use = 4;
% % velocity and height
% figure
% grid on
% [AX1, H11, H12]=plotyy(v_24h(:,1),v_24h(:,2),z_fit(:,1),z_fit(:,2),'stairs','plot');
% set(AX1(1),'xlim',[1 365],'xtick',[1 32 60 91 121 152 182 213 244 274 305 335], ...
% 'xticklabel',{'1/2011','2/2011','3/2011','4/2011','5/2011','6/2011','7/2011' ...
% ,'8/2011','9/2011','10/2011','11/2011','12/2011'},'fontsize',10,'ylim',[0 500],'ytick',[0 100 200 300 400 500],'ycolor','b')
% set(H11,'Color','b','linewidth',1)
% ylabel(AX1(1),{'velocity (my^-^1)'},'fontsize',10)
% title(site,'fontsize',10)
% set(AX1(2),'xlim',[1 365],'xtick',[1 32 60 91 121 152 182 213 244 274 305 335], ...
% 'xticklabel',{'1/2011','2/2011','3/2011','4/2011','5/2011','6/2011','7/2011' ...
% ,'8/2011','9/2011','10/2011','11/2011','12/2011'},'fontsize',10,'ycolor','k')
% ylabel(AX1(2),{'surface';' height (m)'},'fontsize',10)
% set(H12,'Color','k','linewidth',1)