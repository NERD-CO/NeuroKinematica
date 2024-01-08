% filename = fullfile('C:\Users\erinr\Desktop\IO_S1_20230309_test-er-2023-06-29-3d\IO_S1_test-er-2023-06-29\videos\20230309_b2_d0p4_session001_frontCam-0000.mp4');
% filename2 = fullfile('C:\Users\erinr\Desktop\IO_S1_20230309_test-er-2023-06-29-3d\IO_S1_test-er-2023-06-29\videos\20230309_b2_d0p4_session001_frontCam-0000.avi');

%create objects to read and write the video
readerObj = VideoReader(filename);
writerObj = VideoWriter(filename2,'Uncompressed AVI');

%open AVI file for writing
open(writerObj);

%read and write each frame
for k = 1:readerObj.NumFrames
    img = read(readerObj,k);
    writeVideo(writerObj,img);
end
close(writerObj);
