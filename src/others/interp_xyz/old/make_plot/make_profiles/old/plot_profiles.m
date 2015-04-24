% plot cross-section and also point depth profiles

clear all;close all;clc
%% parameter

% load mesh data
load('mesh.dat','-mat','rr','aa','a01','na','nr');

% model data
fin = 'vsh.xyz';

% output
fn_fig = 'vsh.pdf';

%% load model data

data = load(fin);
idx = data(:,1);
v = data(:,5);

vv = nan(size(rr));
for i = 1:numel(idx)
  vv(idx(i)) = v(i);
end

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

pcolor(hax1,yy,xx,vv)
axis equal tight
axis off
shading flat
% colorbar
colormap jet
colormap(flipud(colormap))

%% plot point depth profiles

hold on
scale = 200;
for i = 5:5:na-5
  ai = aa1(1,i);
  vi = vv(:,i);
  vi = vi-5;
  xi = xx(:,i);
  yi = yy(:,i);
  plot(hax1,yi+scale*vi*cos(ai),xi-scale*vi*sin(ai),'k','linewidth',1)
  plot(hax1,yi,xi,'w')
end

%% save figure

formatfig(hf)
print(hf,'-dpdf','-painters',fn_fig)