testplot=true;
nimages=5;
nimplants=100;
maxstreaks=3;

imagepath='/home/enrico/Eran/testimages/';
psffiles=dir(fullfile(imagepath,'*sci_coadd_PSF_1.fits'));
imagefiles=dir(fullfile(imagepath,'*sci_coadd_Image_1.fits'));

% create a table for storing results
detectiontable=table('Size',[nimages*nimplants,5],...
    'VariableTypes',{'string','string','struct','struct','double'},...
    'VariableNames',{'folder','image','implanted','detected','timing'});


for i=1:nimages
    % pick up an image in the image directory and its corresponding PSF
    % (assume that both files are present, so it is sufficient to check for
    %  the existence of one of them)
    ifile=floor(rand*numel(imagefiles))+1;
    AI=AstroImage(fullfile(imagefiles(ifile).folder,imagefiles(ifile).name));
    PSF=AstroImage(fullfile(psffiles(ifile).folder,psffiles(ifile).name));
    AI.PSF=PSF.Image;
    
    % normalize the image
    AI=imProc.background.backVar(AI,'Block',[128 128], 'Method',...
        @imUtil.background.modeVar_LogHist, 'MethodArgs',{{'MinVal',5, 'MaxVal',3000},{}});
    AI=imProc.image.subBackDivideStd(AI);
    
    % implant random streaks
    for j=1:nimplants
        AI1=AI.copy;
        nstreaks=floor(rand*maxstreaks+1);
        streakcoords = [(rand(nstreaks,2)-0.25)*2*size(AI1.Image,1),...
                        (rand(nstreaks,2)-0.25)*2*size(AI1.Image,2)];
        streakstrength = rand(nstreaks,1)*1000; % maybe better log distribution
        for k=1:nstreaks
            AI1.ImageData.Image = imUtil.streaks.addLineToImage(AI1.ImageData.Image, streakcoords(k,:),...
                streakstrength(k), AI.PSFData.Data);
        end
        
        % refilter with PSF (matched filtering?)
        AI1=imProc.image.xcorrWithPSF(AI1);
          
        im=(AI1.Image);
        
        % detect streaks
        tic;
        segs=lsd_scale_mex(im,1/3);
        segs=merge_segments(segs);

        detectiontable.timing(l) = toc;
        
        % store the test results in the table
        l=(i-1)*nimages+j;
        detectiontable.folder(l) = imagefiles(ifile).folder;
        detectiontable.image(l) = imagefiles(ifile).name;
        detectiontable.implanted(l).coords = streakcoords(:,[1 3 2 4]);
        detectiontable.implanted(l).strength = streakstrength;
        detectiontable.detected(l).coords = segs';
        
        % plot for visual feedback
        if testplot
            colormap bone
            imagesc(im,[0,50]); axis ij;
            hold on
            plot(segs([2,4],:),segs([1,3],:),'LineWidth',2)
            hold off
            drawnow
        end
        
    end
end