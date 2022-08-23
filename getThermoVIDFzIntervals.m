% getThermoVIDFzIntervals: 
% Detect freezing periods from thermographic videos (TESTO 885).
% Code used in "Brainstem somatostatin-expressing cells promote analgesia during negative emotional states". Winke, et al. 2022.
% [v0.1, Daniel Jercog 2022, see: github.com/djercog/TestoFreez]
%
% function [fzIntervals,tm,areaChange,areaMouse]=getThermoVIDFzIntervals(videoFullpath,ffmpegBinPath)
% function [fzIntervals,tm,areaChange,areaMouse]=getThermoVIDFzIntervals(videoFullpath,ffmpegBinPath,immobilityThreshold,minFzDur,minNonFzMerge)
%
% Input:
%  videoFullpath: full path to thermographic video (e.g, videoFullpath='E:\Videos\WT_FCA_IR\24-CSP.wmv';)
%  ffmpegBinPath: Path to ffmpeg/bin instalation folder (e.g., ffmpegBinPath='E:\ffmpeg\bin';)
%  immobilityThreshold: Threshold for immobility detection (default immobilityThreshold=0.003;)
%  minFzDur: Minimum immobility period duration to be considered freezing (default minFzDur=0.5;)
%  minNonFzMerge: Minimum non-freezing period to be merged between consecutive freezing periods (default minNonFzMerge=0.1;)
%
% Output:
%  fzIntervals: Freezing intervals detected.
%  tm: Time
%  areaChange: % pixels changed between consecutive video image frames within the ROI defined.
%  areaMouse: mouse area detected in pixels.
%
% Requires: 
% - videoframets (https://www.mathworks.com/matlabcentral/fileexchange/61235-video-frame-time-stamps)
% - ffmpeg (https://ffmpeg.org/download.html)

function [fzIntervals,tm,areaChange,areaMouse]=getThermoVIDFzIntervals(videoFullpath,ffmpegBinPath,immobilityThreshold,minFzDur,minNonFzMerge)

if nargin < 3
    immobilityThreshold = 0.003;
    minFzDur=0.5;
    minNonFzMerge=0.1;
else
    if (nargin~=5)
        error('Wrong number of input arguments.')
    end
end

%check if ffmpeg exists (for vidoframerets0
if exist(ffmpegBinPath)==0
    error('ffmpeg directory not found.')
end
tm=videoframets(ffmpegBinPath,videoFullpath); %get video frame timestamps based on ffmpeg

%set HP roi
vidObj = VideoReader(videoFullpath);
vidFramesAux = read(vidObj,4500);
fAux=figure('color','w');
ax=subplot(1,1,1);hold on;
title('Define HP area:','fontsize',12);
roiHP=roipoly(vidFramesAux);
roiHP=double(roiHP);
roiHP(roiHP==0)=nan;
close(fAux);

%Load entire video
display('Loading video....')
vidFramesAux = read(vidObj);

%Calculate %area change within ROI between binarized downsampled video consecutive frames
display('Processing video....');
clear areaChange areaMouse;
for j=2:length(vidFramesAux)
    %obtain 2 consecutive video frames
    auxF1=roiHP.*double(vidFramesAux(:,:,3,j-1)); %process blue channel only
    auxF2=roiHP.*double(vidFramesAux(:,:,3,j));   %process blue channel only
    thresh=min(auxF1(:))+ graythresh(auxF1)*(max(auxF1(:))-min(auxF1(:))); %Otsu's thresholding
    auxF1(auxF1<=thresh)=0;auxF1(auxF1>thresh)=1;
    thresh=min(auxF2(:))+ graythresh(auxF2)*(max(auxF2(:))-min(auxF2(:))); %Otsu's thresholding
    auxF2(auxF2<=thresh)=0;auxF2(auxF2>thresh)=1;
    %invert scale: 1 on mouse;
    indTo1=find(auxF1==0);indTo0=find(auxF1==1);auxF1(indTo1)=1;auxF1(indTo0)=0;
    indTo1=find(auxF2==0);indTo0=find(auxF2==1);auxF2(indTo1)=1;auxF2(indTo0)=0;
    areaMouse(j-1)=(length(find(auxF1==1))+length(find(auxF2==1)))/2; %avg between the 2 frames
    areaChange(j-1)=length(find(auxF2-auxF1==1 | auxF2-auxF1==-1))/length(find(roiHP==1)); % %change relative to ROI
end

%detecting freezing
tm=tm(2:end);
fzAux=zeros(length(areaChange),1);
fzAux(find(areaChange<immobilityThreshold))=1;
endFzPer=find(diff(fzAux)==-1);
begFzPer=find(diff(fzAux)==1)+1;
if ~isempty(begFzPer) && ~isempty(endFzPer)
    if endFzPer(1)<begFzPer(1)
        begFzPer=[1;begFzPer];
    end
    if begFzPer(end)>endFzPer(end)
        endFzPer=[endFzPer;length(fzAux)];
    end
    fzOns=tm(begFzPer);
    fzOffs=tm(endFzPer);
else
    fzOns=[];
    fzOffs=[];
end
%remove short nFZ periods if neighboring periods are fz periods
indAux=find(fzOns(2:end)-fzOffs(1:end-1)<=minNonFzMerge);
toRemove=[];
for s=1:length(indAux)
        %if (fzOffs(indAux(s)-1)-fzOns(indAux(s)-1)>minFzDur && fzOffs(1+indAux(s))-fzOns(1+indAux(s))>minFzDur)
        if (fzOffs(1+indAux(s))-fzOns(1+indAux(s))>minNonFzMerge)
            toRemove=[toRemove,indAux(s)];
        end
end
fzOffs(toRemove)=[];
fzOns(toRemove+1)=[];
%remove shor fz
indAux=find(fzOffs-fzOns<=minFzDur);
fzOns(indAux)=[];
fzOffs(indAux)=[];
fzIntervals=[fzOns,fzOffs];

verifyDet = input('Verify detection? Y/N [N]:','s');
if strcmp(verifyDet,'Y')
    % Verify detection
    f = figure('color','w');
    ax(1)=subplot(5,1,5);hold on;
    plot(tm,areaChange,'k');
    axis tight;ylim([-0.005,0.02]);
    plot(fzIntervals,[0.01,0.01],'r','linewidth',2);
    ylabel('% Area Change');xlabel('Time (s)');
    set(gca,'tickdir','out');
    ax(2)=subplot(5,1,1:4,'fontsize',5);hold on;
    title('Verify detection...','fontsize',14);
    axis ij;
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
    for i=1:size(vidFramesAux,4)-1
        if ishandle(f)
            im=imagesc(vidFramesAux(:,:,:,i));
            p=plot(ax(1),[tm(i),tm(i)],[-0.005,0.02],'m','linewidth',1);
            pause(0.02);
            delete(im);
            delete(p);
        end
    end
    close(f);
end