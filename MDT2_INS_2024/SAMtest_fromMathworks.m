I = imread("peppers.png");
imageSize = size(I);

%%

I = 255.*rescale(I);

%%

model = segmentAnythingModel;

%%

embeddings = extractEmbeddings(model,I);

%%

pointPrompt = [453 283; 496 288; 504 300];
backgroundPoints = [308 272; 348 176];

%%

[masks,scores,maskLogits] =  segmentObjectsFromEmbeddings(model,embeddings,imageSize, ...
    ForegroundPoints=pointPrompt,BackgroundPoints=backgroundPoints);

%%

[masks2,scores2,maskLogits2] =  segmentObjectsFromEmbeddings(model,embeddings,imageSize, ...
    ForegroundPoints=pointPrompt,BackgroundPoints=backgroundPoints,MaskLogits=maskLogits);

%%
figure;
imMask1 = insertObjectMask(I,masks);
imshow(imMask1)

%%

figure;
imMask2 = insertObjectMask(I,masks2);
imshow(imMask2)