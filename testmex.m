% nx=1700; ny=1700;
% AI=AstroImage();
% AI.Image=2000*imUtil.art.createSegments([nx,ny],...
%     [322,2033;98,0],[540,11;145,211],'width',1.5) + ...
%     1000*rand(nx,ny);
AI=AstroImage('/home/enrico/Eran/testimages/LAST.01.02.01_20251203.030427.430_clear_1423_000_001_014_sci_coadd_Image_1.fits');

PSF=AstroImage('/home/enrico/Eran/testimages/LAST.01.02.01_20251203.030427.430_clear_1423_000_001_014_sci_coadd_PSF_1.fits');
AI.PSF=PSF.Image;

AI=imProc.background.backVar(AI,'Block',[128 128], 'Method',...
    @imUtil.background.modeVar_LogHist, 'MethodArgs',{{'MinVal',5, 'MaxVal',3000},{}});
AI=imProc.image.subBackDivideStd(AI);

AI.ImageData.Image = imUtil.streaks.addLineToImage(AI.ImageData.Image, [1345,678,998,109],...
                                   30, AI.PSFData.Data);
AI.ImageData.Image = imUtil.streaks.addLineToImage(AI.ImageData.Image, [1845,478,1098,509],...
                                   20, AI.PSFData.Data);

AI=imProc.image.xcorrWithPSF(AI);


im=double(AI.Image);
%im(im>100)=100;

%im=double(rgb2gray(imread('../Gary_LSD/code/undistortedImage/1.png')));
%im =double(AI.Image);
%im=cumsum(AI.Image); im=im-repmat(sum(im,2)/ny,1,ny);
%im=cumsum(AI.Image-mean(AI.Image(:)),2);
%im = double(imread('cameraman.tif'));
[Hxx,Hxy,Hyy,lambda1,lambda2]=hessian(im);
%im=log(-lambda2.*(lambda2<-1))*9000;

tic;
segs=lsd_scale_mex(im,1/3);
toc

colormap bone
imagesc(im,[0,50]); axis ij;
%imagesc(log(-lambda2).*(lambda2<0)); axis ij; colorbar
hold on
plot(segs([2,4],:),segs([1,3],:),'LineWidth',2)
hold off

