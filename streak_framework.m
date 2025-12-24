testplot=true;
nimages=5;
nimplants=100;

% create a table for storing results
detectiontable=table;

for i=1:nimages
    % pick up an image in the image directory and its corresponding PSF
    % (assume that both files are present, so it is sufficient to check for
    %  the existence of one of them)
    imagepath='/home/enrico/Eran/testimages/';
    
    AI=AstroImage('/home/enrico/Eran/testimages/LAST.01.02.01_20251203.030427.430_clear_1423_000_001_014_sci_coadd_Image_1.fits');
    
    PSF=AstroImage('/home/enrico/Eran/testimages/LAST.01.02.01_20251203.030427.430_clear_1423_000_001_014_sci_coadd_PSF_1.fits');
    AI.PSF=PSF.Image;
    
    % normalize the image
    
    AI=imProc.background.backVar(AI,'Block',[128 128], 'Method',...
        @imUtil.background.modeVar_LogHist, 'MethodArgs',{{'MinVal',5, 'MaxVal',3000},{}});
    AI=imProc.image.subBackDivideStd(AI);
    
    % implant random streaks
    for j=1:nimplants
        nstreaks=floor(rand*4);
        streakcoords = rand(nstreaks,4)*2*size(AI.Image,1);
        streakstrength = rand*1000; % maybe better log distribution
        for k=1:nstreaks
            AI.ImageData.Image = imUtil.streaks.addLineToImage(AI.ImageData.Image, streakcoords(k,:),...
                streakstrength(k), AI.PSFData.Data);
        end
        
        AI=imProc.image.xcorrWithPSF(AI);
          
        im=(AI.Image);
        
        % detect streaks
        tic;
        segs=lsd_scale_mex(im,1/3);
        toc
        
        segs=merge_segments(segs);
        
        % store the test results in the table
        
        % plot for visual feedback
        if testplot
            colormap bone
            imagesc(im,[0,50]); axis ij;
            hold on
            plot(segs([2,4],:),segs([1,3],:),'LineWidth',2)
            hold off
        end
        
    end
end