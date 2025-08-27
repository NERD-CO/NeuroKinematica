function compressVideosJAT(filedir)
if ~exist('filedir','var')
    filedir = uigetdir('Choose data directory', 'Choose data directory');
end
vidList = getVidList(filedir);
for v = 1:numel(vidList)
    vH = VideoReader(vidList{v});
    destPath = cat(2,vidList{v}(1:end-4),'.mp4');
    wH = VideoWriter(destPath,'MPEG-4');
    wH.Quality = 100;
    wH.open()
    while vH.hasFrame()
        I = vH.readFrame();
        writeVideo(wH,I)
    end
    wH.close()
    [~,vidName] = fileparts(destPath);
    try
        vChk = VideoReader(destPath);
        if vChk.NumberOfFrames == vH.NumberOfFrames
            delete(vH)
            delete(vidList{v})
            disp(cat(2,'Successfully compressed ',vidName))
        end
    catch ME
        getReport(ME)
        disp(cat(2,'Compression error: ',vidName))
    end
end
end

function vidList = getVidList(dirIn)
[files, folders] = rdir(dirIn);
vidList = cell(0,1);
for y = 1:numel(files)
    if contains(files{y},'.avi')
        vidList = cat(1,vidList,files(y));
    end
end
% for x = 1:numel(folders)
%     %Check to make sure that this part of 'a' is an actual file or
%     %subdirectory.
%     if ~strcmp(folders{x},'.')&&~strcmp(folders{x},'..')
%         %Ooooooh, recursive!  Fancy!
%         vidList = cat(1,vidList,getVidList(files{x}));
%     end
% end
end

function [files,folders,size] = rdir(absDir)
%Returns paths of all subfolders, files, subfolder files in a given
%directory, and size of the directory
%Input - (absDir) Absolute path to the directory of interest, ending in \
%Output - [files,folders,size] files: all files located in all absDir and
%subdirectories, folders: all folders located in all absDir and
%subdirectories, size: size in bytes of the absolute directory
fileList=dir(absDir);
folders={};
files={};
tempFiles={};
tempFolders={};
size=0;
if ~strcmp('\',absDir(end))
    absDir=[absDir '\'];
end
for i=1:length(fileList)
    if ~strcmp(fileList(i).name,'..') && ~strcmp(fileList(i).name,'.') && ~strncmp(fileList(i).name,'$',1)
        if isfolder([absDir fileList(i).name])
            folders{end+1,1}=[absDir fileList(i).name '\'];
            if isempty(strfind([absDir fileList(i).name '\'],'.nff'))
                [tempFiles, tempFolders,tempSize]=rdir([absDir fileList(i).name '\']);
            end
            if ~isempty(tempFiles)
                files=[files; tempFiles];
                size=size+tempSize;
            end
            if ~isempty(tempFolders)
                folders=[folders; tempFolders];
            end
        else
            files{end+1,1}=[absDir fileList(i).name];
            size=size+fileList(i).bytes;
        end
    end
end
end