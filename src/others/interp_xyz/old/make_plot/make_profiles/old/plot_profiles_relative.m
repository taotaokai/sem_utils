% plot cross-section (relative) and also point depth profiles

clear all;close all;clc
%% parameter

% load mesh data
load('mesh.dat','-mat','rr','aa','a01','na','nr','r');

% model data
fin = 'vsh.xyz';

% Reference Earth Model (STW105)
%radius[m] density[kg/m^3] vpv[m/s] vsv[m/s] Q_kappa Q_miu vph[m/s] vsh[m/s] eta[m/s]
fn_REF = 'STW105.txt';

% output
fn_fig = 'vsh_stw105.pdf';

%% load Reference Earth Model

fid = fopen(fn_REF,'r');
em = textscan(fid,'%f %f %f %f %f %f %f %f %f %*[^\n]','commentstyle','#');
fclose(fid);

em_r = em{1}/1000;
em_vsh = em{8}/1000;

%% load model data

data = load(fin);
idx = data(:,1);
v = data(:,5);

vv = nan(size(rr));
for i = 1:numel(idx)
  vv(idx(i)) = v(i);
end

%% get reference model at mesh grids

% vv_ref = nan(size(rr));
vv_ref = interp1(em_r,em_vsh,rr);

%% define figure

hf = figure;
% hf = figure('visible','off');

% figure size
set(hf,'paperunits','centimeters','units','centimeters',...
    'PaperPositionMode','manual','papertype','a4',...
    'paperorientation','landscape',...
    'position',[0 0 29.7 21],'paperposition',[0 0 29.7 21])

% paper margins
West = 1.5;
South = 2;
Width = 29.7-2;
Height = 21-5;
% Width = 21-4;
% Height = 29.7-5;

% define axes positions
pos_ax1 = [West,South,Width,Height];
hax1 = axes('Units','centimeters','position',pos_ax1);

%% plot cross-section

aa1 = aa-a01/2;
[xx,yy] = pol2cart(aa1,rr);

dlnvv = (vv-vv_ref)./vv_ref;

pcolor(hax1,yy,xx,100*dlnvv)
axis equal tight
shading flat
% colorbar
colormap jet
colormap(flipud(colormap))
caxis([-6 6])

%% plot point depth profiles

hold on
scale = 200;
for i = 5:5:na-5
  ai = aa1(1,i);
  vi = vv(:,i);
%   vi = vi-mean(vi(~isnan(vi)));
  vi = vi-5;
  xi = xx(:,i);
  yi = yy(:,i);
  plot(hax1,yi+scale*vi*cos(ai),xi-scale*vi*sin(ai),'k','linewidth',1)
  plot(hax1,yi,xi,'w')
end

%% save figure

formatfig(hf)
print(hf,'-dpdf','-painters',fn_fig)