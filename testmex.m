nx=1700; ny=1700;
AI=AstroImage();
AI.Image=2000*imUtil.art.createSegments([nx,ny],...
    [322,2033;98,0],[540,11;145,211],'width',1.5) + ...
    1000*rand(nx,ny);
imProc.background.background(AI);

im=double(rgb2gray(imread('../Gary_LSD/code/undistortedImage/1.png')));
%im =AI.Image;
%im=cumsum(AI.Image); im=im-repmat(sum(im,2)/ny,1,ny);
%im = double(imread('cameraman.tif'));

tic;
segs=lsd_mex(im);
toc

colormap bone
imagesc(im); axis ij;
hold on
plot(segs([2,4],:),segs([1,3],:),'r')
hold off

