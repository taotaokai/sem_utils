function SACST_vespagram(sacdir,saclst,t0,t1,fs,...
                         sachd_rayp,nthroot,str_title,fignm)
% SACST_vespagram

% clear all;close all;clc
%% parameters

%station info
% stnm = 'YP.NE86';

% data
% sacdir = 'RF/';
% saclst = [sacdir,'rayp.bin.lst'];

%time window
% t0 = -10;
% t1 = 30;
% fs = 40;
dt = 1/fs;
t = t0:dt:t1;
nt = length(t);

%sachd: ray parameter 
% sachd_rayp = 'user0';

%scan time window (moveout curve: t=t0+k*p^2)
ksc = -600:400;

%stacking
% nthroot = 2;

%output figure
% str_title = sprintf('%s: vepsagram (%d-th root stacking)',stnm,nthroot);
% fignm = sprintf('%s.vepsagram.pdf',stnm);

%% load sac

sacst = SACST_fread('list',saclst,'prefix',sacdir);
sacst = SACST_interp(sacst,t,'ref','0');
nsac = length(sacst);
comm = sprintf('[sacst.%s]',sachd_rayp);
rayp = eval(comm);

%% vespagram

nk = length(ksc);
vesp = zeros(nk,nt);
x = [sacst.data];
beam = zeros(nt,nsac);
for ik = 1:nk
    tshift = ksc(ik)*rayp.^2;
    beam(:) = 0;
    for isac = 1:nsac
        xi = interp1(t,x(:,isac),t+tshift(isac));
        xi(isnan(xi)) = 0;
        beam(:,isac) = xi;
    end
    % n-th root stacking
    beam = mean(sign(beam).*abs(beam).^(1/nthroot),2);
    beam = sign(beam).*abs(beam).^nthroot;
    vesp(ik,:) = beam/nsac;
end

%% Define plot

hf = figure;
% hf = figure('visible','off');

% figure size
set(hf,'paperunits','centimeters','units','centimeters',...
    'PaperPositionMode','manual','papertype','a4',...
    'paperorientation','portrait',...
    'position',[0 0 21 29.7],'paperposition',[0 0 21 29.7])

% paper margins
West = 2;
South = 2;
Width = 21-4;
% Height = 29.7-4;

% define axes positions
pos_ax = [West,South,Width,Width/2];
hax = axes('Units','centimeters','position',pos_ax);
%colorbar
cbarwidth = 0.2;
cbarheight = Width/4;
pos_cbar = [West+Width+0.2,South+Width/8,cbarwidth,cbarheight];

%% vespagram

set(hf,'currentaxes',hax)
vesp = vesp/max(vesp(:));
data = abs(vesp);
imagesc(t,ksc,data); set(hax,'ydir','normal')
caxis([0 0.3])
hold on 
plot([t0 t1],[0 0],'w','linewidth',2);
% plot([0 0],[min(ksc) max(ksc)],'w','linewidth',2)
xlabel('Time (sec)')
ylabel('Moveout coefficient (km^2/s)')
title(str_title)

%colorbar
hcbar = colorbar('location','eastoutside','Units','centimeters',...
        'position',pos_cbar);
set(hcbar,'Units','centimeters',...
    'ticklength',[cbarwidth/cbarheight,1])
ylabel(hcbar,'normalized stack amplitude')
set(hax,'Units','centimeters','position',pos_ax);

formatfig(hf)
print(hf,'-dpdf','-painters',fignm)

end