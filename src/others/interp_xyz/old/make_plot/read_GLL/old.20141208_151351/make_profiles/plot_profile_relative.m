function plot_profile_relative(profilename, modelname, scalefactor)
% plot cross-section (relative) and also point depth profiles
% clear all;close all;clc
%% parameter

% load mesh data
load('mesh/mesh.dat','-mat','rr','aa','a01','na','nr');

% Reference Earth Model (STW105)
%radius[m] density[kg/m^3] vpv[m/s] vsv[m/s] Q_kappa Q_miu vph[m/s] vsh[m/s] eta[m/s]
fn_REF = 'REFmodel/STW105.txt';

% model data
fin = sprintf('model/%s_%s.xyz',profilename,modelname);

% output
fn_fig = sprintf('figure/%s_%s_STW105.pdf',profilename,modelname);

% scale factor for point depth profile
if nargin < 2
  scalefactor = 100;
end

% color range
clim = [-6,6];

%% load Reference Earth Model
fid = fopen(fn_REF,'r');
em = textscan(fid,'%f %f %f %f %f %f %f %f %f %*[^\n]','commentstyle','#');
fclose(fid);

em_r = em{1}/1000;
%radius[m] density[kg/m^3] vpv[m/s] vsv[m/s] 
%Q_kappa Q_miu vph[m/s] vsh[m/s] eta[m/s]
switch modelname
  case 'rho'
    em_v = em{2}/1000;
  case 'vpv'
    em_v = em{3}/1000;
  case 'vsv'
    em_v = em{4}/1000;
  case 'vph'
    em_v = em{7}/1000;
  case 'vsh'
    em_v = em{8}/1000;
  case 'eta'
    em_v = em{9};
end

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
vv_ref = interp1(em_r,em_v,rr);

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

%colorbar
cbarheight = 0.2;
cbarwidth = Width/4;
pos_cbar = [West+Width/2-cbarwidth/2,South+3,cbarwidth,cbarheight];

%% plot cross-section

aa1 = aa-a01/2;
[xx,yy] = pol2cart(aa1,rr);

dlnvv = (vv-vv_ref)./vv_ref;

pcolor(hax1,yy,xx,100*dlnvv)
axis equal tight
shading flat
caxis(clim)

% set colorbar
hcbar = colorbar('location','southoutside','Units','centimeters',...
        'position',pos_cbar);
set(hcbar,'Units','centimeters','ticklength',[cbarheight/cbarwidth,1])
xlabel(hcbar,modelname)
colormap jet
colormap(flipud(colormap))

set(hax1,'Units','centimeters','position',pos_ax1);

%% plot point depth profiles
hold on

v_mean = mean(vv,2);
v_mean = mean(v_mean(~isnan(v_mean)));
for i = 5:5:na-5
  ai = aa1(1,i);
  vi = vv(:,i);
  vi = vi-v_mean;
%   vi = vi-5;
  xi = xx(:,i);
  yi = yy(:,i);
  plot(hax1,yi+scalefactor*vi*cos(ai),...
       xi-scalefactor*vi*sin(ai),'k','linewidth',1)
  plot(hax1,yi,xi,'w')
end

%% save figure

formatfig(hf)
print(hf,'-dpdf','-painters',fn_fig)

end