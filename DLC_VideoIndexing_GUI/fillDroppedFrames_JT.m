function fillDroppedFrames_JT(path2videos)
% Get video list in path2videos directory
vidListMP4 = dir(fullfile(path2videos,'*.mp4'));
vidListMP4 = {vidListMP4(:).name};
vidListAVI = dir(fullfile(path2videos,'*.avi'));
vidListAVI = {vidListAVI(:).name};
vidList = cat(1,vidListMP4(:),vidListAVI(:));
% Pre-allocate
tDrop = cell(numel(vidList),1);
framesWereDropped = false(numel(vidList),1);
maxfrms = zeros(numel(vidList),1);
% Load the timestamps and find dropped frames
for vndx = 1:numel(vidList)
    timePath = fullfile(path2videos,cat(2,vidList{vndx}(1:end-4),'_timestamps.txt'));
    if isfile(timePath)
        timestamps = importdata(timePath);
        timeInts = round(timestamps/median(timestamps(5:100)));
        if max(timeInts(:)) > 1
            framesWereDropped(vndx) = true;
        end
        tDrop{vndx} = cumsum(timeInts);
        maxfrms(vndx) = tDrop{vndx}(end);
    end
end
for vndx = 1:numel(vidList)
    % Define read and write paths
    vidPath = fullfile(path2videos,vidList{vndx});
    destPath = fullfile(path2videos,cat(2,vidList{vndx}(end-2:end),'_dropsFilled.mp4'));
    % If frames were dropped or video needs compressing or videos are different lengths:
    if framesWereDropped(vndx) || contains(vidPath,'.avi') || diff(maxfrms) > 0
        % Initialize video read/write objects
        hV = VideoReader(vidPath);
        hW = VideoWriter(destPath,'MPEG-4');
        hW.Quality = 100;
        open(hW)
       
       
        frmReadCt = 0;
        frmWriteCt = 0;
        while hV.hasFrame()
            frm = readFrame(hV);
            frmReadCt = frmReadCt+1;
            for fml = 1:tDrop{vndx}(frmReadCt)
                writeVideo(hW,frm);
                frmWriteCt = frmWriteCt+1;
            end
           
        end
        for endFill = 1:(max(maxfrms)-frmWriteCt)
            writeVideo(hW,frm);
        end
        close(hW)
    else
        copyfile(vidPath,destPath)
    end
end