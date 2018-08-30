clear all;
close all;

y=load('image_sc.txt');

figure(1)
surf(y);
shading interp;
colormap gray;
view([0 90]);
axis off;