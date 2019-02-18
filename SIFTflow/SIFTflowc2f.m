% function for coarse to fine SIFT flow matching
% 
% Usage:
%
% [vx,vy,energylist]=SIFTflowc2f(im1,im2,SIFTflowpara,isdisplay);
%
% Input arguments:
%     im1 an im2: SIFT images to match (from im1 to im2). The two images are not required
%                 to have the same dimension, but it is not likely to have meaningful 
%                 matches if the two images are too different in size
%
%     SIFTflowpara: a structure parameter for SIFT flow matching. Here are the fields
%           alpha:   (0.01) the coefficient of the truncated L1-norm regularization of the flow discontinuity
%           d:       (alpha*20) the threshold of the truncated L1-norm
%           gamma:   (0.001) the coefficient of the regularization on flow magnitude (L1 norm)
%           nlevels: (4) the number of levels of the Gaussian pyramid
%           topwsize: (10) the size of the matching window at the top level
%           nTopIterations: (100) the number of BP iterations at the top level
%           wsize:   (3) the size of the matching window at lower levels
%           nIterations: (40) the number of BP iterations at lower levels
%
%     isdisplay: (false) a boolean variable indicating whether information should be displayed
%
% Ce Liu
% celiu@mit.edu
% CSAIL MIT
% Jan, 2009
%              
% All Rights Reserved

function [vx,vy,energylist]=SIFTflowc2f(im1,im2,SIFTflowpara,isdisplay)

if isfield(SIFTflowpara,'alpha')
    alpha=SIFTflowpara.alpha;
else
    alpha=0.01;
end

if isfield(SIFTflowpara,'d')
    d=SIFTflowpara.d;
else
    d=alpha*20;
end

if isfield(SIFTflowpara,'gamma')
    gamma=SIFTflowpara.gamma;
else
    gamma=0.001;
end

if isfield(SIFTflowpara,'nlevels')
    nlevels=SIFTflowpara.nlevels;
else
    nlevels=4;
end

if isfield(SIFTflowpara,'wsize')
    wsize=SIFTflowpara.wsize;
else
    wsize=3;
end

if isfield(SIFTflowpara,'topwsize')
    topwsize=SIFTflowpara.topwsize;
else
    topwsize=10;
end

if isfield(SIFTflowpara,'nIterations')
    nIterations=SIFTflowpara.nIterations;
else
    nIterations=40;
end

if isfield(SIFTflowpara,'nTopIterations')
    nTopIterations=SIFTflowpara.nTopIterations;
else
    nTopIterations=100;
end

if exist('isdisplay','var')~=1
    isdisplay=false;
end
if ~isfloat(im1)
    im1=im2double(im1);
end
if ~isfloat(im2)
    im2=im2double(im2);
end

% build the Gaussian pyramid for the SIFT images
pyrd(1).im1=im1;
pyrd(1).im2=im2;

for i=2:nlevels
    pyrd(i).im1=imresize(imfilter(pyrd(i-1).im1,fspecial('gaussian',5,0.67),'same','replicate'),0.5,'bicubic');
    pyrd(i).im2=imresize(imfilter(pyrd(i-1).im2,fspecial('gaussian',5,0.67),'same','replicate'),0.5,'bicubic');
end

for i=1:nlevels
    [height,width,nchannels]=size(pyrd(i).im1);
    [height2,width2,nchannels]=size(pyrd(i).im2);
    [xx,yy]=meshgrid(1:width,1:height);    
    pyrd(i).xx=round((xx-1)*(width2-1)/(width-1)+1-xx);
    pyrd(i).yy=round((yy-1)*(height2-1)/(height-1)+1-yy);
end

nIterationArray=round(linspace(nIterations,nIterations*0.6,nlevels));

for i=nlevels:-1:1
    if isdisplay
        fprintf('Level: %d...',i);
    end
    [height,width,nchannels]=size(pyrd(i).im1);
    [height2,width2,nchannels]=size(pyrd(i).im2);
    [xx,yy]=meshgrid(1:width,1:height);

    if i==nlevels % top level
        vx=pyrd(i).xx;
        vy=pyrd(i).yy;
        
        winSizeX=ones(height,width)*topwsize;
        winSizeY=ones(height,width)*topwsize;
    else % lower levels
        vx=round(pyrd(i).xx+imresize(vx-pyrd(i+1).xx,[height,width],'bicubic')*2);
        vy=round(pyrd(i).yy+imresize(vy-pyrd(i+1).yy,[height,width],'bicubic')*2);
        
        winSizeX=ones(height,width)*wsize;
        winSizeY=ones(height,width)*wsize;
    end
    if nchannels<=3
        Im1=im2feature(pyrd(i).im1);
        Im2=im2feature(pyrd(i).im2);
    else
        Im1=pyrd(i).im1;
        Im2=pyrd(i).im2;
    end
    
    if i==nlevels
        [flow,foo]=mexDiscreteFlow(Im1,Im2,[alpha,d,gamma*2^(i-1),nTopIterations,2,topwsize],vx,vy,winSizeX,winSizeY);
    else
        [flow,foo]=mexDiscreteFlow(Im1,Im2,[alpha,d,gamma*2^(i-1),nIterationArray(i),nlevels-i,wsize],vx,vy,winSizeX,winSizeY);
    end
    energylist(i).data=foo;
    vx=flow(:,:,1);
    vy=flow(:,:,2);
    if isdisplay
        fprintf('done!\n');
    end
end




