%% SETUP FOR PLOTS

% Declare this top header once.

%Sets the units of your root object (screen) to centimeters
set(0,'units','centimeters');
%Obtains this inch information
fig_setup = struct;

fig_setup.CM_SS = get(0,'screensize');
fig_setup.width = fig_setup.CM_SS(3);
fig_setup.height= fig_setup.CM_SS(4);
fig_setup.fig_wd = 9.8; %width of figure in cm
fig_setup.fig_wd_wide = 21;
fig_setup.fig_hgt= 5; %height of figure in cm
fig_setup.fig_hgt_wide= 7.5;
fig_setup.fntsize = 8; %fontsize in px
 
fig_setup.export_figures = 1;
% 'vector' for pdf
% 'image'  for png
fig_setup.img_format = 'vector';

fig_setup.color(1) = "#0047AB"; %cobalt blue
fig_setup.color(2) = "#FF5733"; %tomato
fig_setup.color(3) = "#13BF1D"; %dark pastel green
fig_setup.color(4) = "#A41495"; %purple
fig_setup.color(5) = "#E9E16D"; %dark yellow
fig_setup.color(6) = "#838383"; %dark grey
fig_setup.color(7) = "#87A5CF"; %blue-light

if strcmp(fig_setup.img_format,'vector')
    fig_setup.img_ext = ".pdf";
else
    fig_setup.img_ext = ".png";
end