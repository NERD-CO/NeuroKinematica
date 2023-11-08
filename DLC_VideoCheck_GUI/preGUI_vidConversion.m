path = 'C:\Users\erinr\Desktop\DLC';
video = "20230912_idea08_session022_rightCam-0000DLC_resnet50_Clin_2023-09-12_LSTN_v3Oct3shuffle1_100000_labeled.mp4";

tab_vidObj = VideoReader(video);
            tab_vid = struct('cdata',zeros(tab_vidObj.Height,tab_vidObj.Width,3,'uint8'),'colormap',[]);

            frami = 1;
            while hasFrame(tab_vidObj)
                tab_vid(frami).cdata = readFrame(tab_vidObj);
                frami = frami+1;
            end
            app.LEFTdlcMainData = struct2table(tab_vid);
