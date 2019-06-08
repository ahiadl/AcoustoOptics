close all;
clear all;
clc;

fullImage = imread('./GUI/arrows.png');

left  = fullImage(1:116, 1:116,  :);
right = fullImage(1:116, 162:277,:);
up    = fullImage(1:116, 324:439,:);
down  = fullImage(1:116, 485:end,:);

IconStack(:,:,1:3)   = left; 
IconStack(:,:,4:6)   = right;
IconStack(:,:,7:9)   = up;
IconStack(:,:,10:12) = down;

names= {'Left'; 'Right'; 'Up'; 'Down'};

for i=1:4
    im = IconStack(:,:,(i-1)*3+1:i*3);
    [row, col] = find(rgb2gray(im)==0);
    n = length(row);
    idxs = [row, col,   ones(n,1);
            row, col, 2*ones(n,1);
            row, col, 3*ones(n,1)];  
    subs = sub2ind(size(left), idxs(:,1), idxs(:,2),idxs(:,3));
    im (subs) = floor(256*0.96);
    imwrite(im,['./GUI/', names{i},'.png']);
end


