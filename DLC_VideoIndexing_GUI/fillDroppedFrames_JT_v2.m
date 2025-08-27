function fillDroppedFrames_JT_v2(path2videos, quality)
    if nargin < 2
        quality = 100;
    end
    
    % Initialize parameters and pre-allocate variables
    [vidList, numVideos] = getVideoList(path2videos);
    [tDrop, framesWereDropped, maxfrms] = preAllocateMemory(numVideos);

    % Load the timestamps and find dropped frames
    [tDrop, framesWereDropped, maxfrms] = findDroppedFrames(path2videos, vidList, numVideos, tDrop, framesWereDropped, maxfrms);
    
    % Fill in the dropped frames and write new videos
    writeNewVideos(path2videos, vidList, numVideos, tDrop, framesWereDropped, maxfrms, quality);
end

% Get video list in path2videos directory
function [vidList, numVideos] = getVideoList(path2videos)
    vidListMP4 = dir(fullfile(path2videos,'*.mp4'));
    vidListAVI = dir(fullfile(path2videos,'*.avi'));
    vidList = cat(1, {vidListMP4(:).name}, {vidListAVI(:).name});
    numVideos = numel(vidList);
end

% Pre-allocate memory
function [tDrop, framesWereDropped, maxfrms] = preAllocateMemory(numVideos)
    tDrop = cell(numVideos, 1);
    framesWereDropped = false(numVideos, 1);
    maxfrms = zeros(numVideos, 1);
end

% Find the timestamps and identify dropped frames
function [tDrop, framesWereDropped, maxfrms] = findDroppedFrames(path2videos, vidList, numVideos, tDrop, framesWereDropped, maxfrms)
    for vndx = 1:numVideos
        timePath = fullfile(path2videos, strcat(vidList{vndx}(1:end-4),'_timestamps.txt'));
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
end

% Write the new videos, filling in any dropped frames
function writeNewVideos(path2videos, vidList, numVideos, tDrop, framesWereDropped, maxfrms, quality)
    for vndx = 1:numVideos
        % Define read and write paths
        vidPath = fullfile(path2videos, vidList{vndx});
        destPath = fullfile(path2videos, strcat(vidList{vndx}(1:end-4),'_dropsFilled.mp4'));
        
        % Check conditions for frame fill or copy
        if framesWereDropped(vndx) || contains(vidPath,'.avi') || any(diff(maxfrms))
            % Initialize video read/write objects
            hV = VideoReader(vidPath);
            hW = VideoWriter(destPath,'MPEG-4');
            hW.Quality = quality;
            open(hW);
            
            frmReadCt = 0;
            frmWriteCt = 0;
            while hV.hasFrame()
                frm = readFrame(hV);
                frmReadCt = frmReadCt + 1;
                for fml = 1:tDrop{vndx}(frmReadCt)
                    writeVideo(hW,frm);
                    frmWriteCt = frmWriteCt + 1;
                end
            end
            
            % Fill in extra frames to make all videos the same length
            for endFill = 1:(max(maxfrms) - frmWriteCt)
                writeVideo(hW,frm);
            end
            
            close(hW);
        else
            copyfile(vidPath, destPath);
        end
    end
end
