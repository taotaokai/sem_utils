% compare tripllicated SH wave with specfem synthetics
clear all;close all;clc

%% parameter

% data
datdir = 'bp/';
list = [datdir,'BXT.lst'];
cmpnm = {'.BXT'};

% travel-time info
tt_ab = load('ttimes/S_660ab');
tt_bc = load('ttimes/S_660bc');
tt_cd = load('ttimes/S_660cd');

% plot time window
t0 = 50;
t1 = 200;
dt = 0.01;
reduce_rayp = 16; % s/deg
ts = t0:dt:t1;

% trace amplification coeff.
tramp = 0.5;

% shift synthetics
t_shift = 0; % positive: shift to positive time 

% polarity
polarity = -1; % due to reversed y axis direction

% output
eventid = '20030727_0625';
fignm = [eventid,'_BXT.pdf'];

%% read sac

sacst = SACST_fread('list',list, 'prefix',datdir, 'suffix',cmpnm);

nsta = size(sacst,1);
gcarc = [sacst(:,1).gcarc];

%% plot

hf = figure;
% hf = figure('visible','off');

% figure size
set(hf,'paperunits','centimeters','units','centimeters',...
    'PaperPositionMode','manual','papertype','a4',...
    'paperorientation','portrait',...
    'position',[0 0 21 29.7],'paperposition',[0 0 21 29.7])

% paper margins
West = 1.5;
South = 2;
Width = 21-4;
Height = 29.7-5;

% define axes positions
% seismic traces
pos_ax1 = [West,South,Width,Height];
hax1 = axes('Units','centimeters','position',pos_ax1);

for i = 1:nsta
    t_reduce = gcarc(i)*reduce_rayp;
    sacsti = SACST_interp(sacst(i,:),t_reduce+ts,'hdr','0');
    syn = sacsti(1).data;
    % normalize amplitude
    syn = syn/max(abs(syn));
    %
    plot(ts+t_shift,polarity*tramp*syn+gcarc(i),'r'); hold on
    % plot station name
    staid = sprintf(' %s(%.1f)',deblank(sacsti(1).kstnm),sacsti(1).az);
    text(t1,gcarc(i),staid,'verticalalignment','top');
end

% plot travel-time 
plot(tt_ab(:,2)-reduce_rayp*tt_ab(:,1), tt_ab(:,1),'b');
plot(tt_bc(:,2)-reduce_rayp*tt_bc(:,1), tt_bc(:,1),'b');
plot(tt_cd(:,2)-reduce_rayp*tt_cd(:,1), tt_cd(:,1),'b');

xlim([t0,t1]);
ylim([min(gcarc)-1, max(gcarc)+1])
set(hax1,'ydir','reverse')

str_xlabel = sprintf('t-dist*%d (s)',reduce_rayp);
xlabel(str_xlabel);
str_ylabel = sprintf('dist (degree)');
ylabel(str_ylabel)
str_title = sprintf('EventID %s SH',eventid);
title(str_title,'interpreter','none')
legend('syn')
legend('boxoff')

%% save figure

formatfig(hf)
print(hf,'-dpdf','-painters',fignm)

% end