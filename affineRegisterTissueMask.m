function newTissueMask = affineRegisterTissueMask(targetModel,matchModel,scalpData)
%This function scales the voxelised tissue mask to the size of the scaled head model.
%This function relies on the DOT-HUB toolbox.

if isequal(scalpData,'landmarks')
    dataForTransTarget = targetModel.landmarks;
    dataForTransMatch = matchModel.landmarks;
elseif isequal(scalpData,'10_20')
    idx_10_20 = [1,6,13,18,22,29,33,35,39,46,50,53,57,64,68,73,80,87,94];
    dataForTransTarget = targetModel.refpts_10_5(idx_10_20,:);
    dataForTransMatch = matchModel.refpts_10_5(idx_10_20,:);
elseif isequal(scalpData,'10_10')
    idx_10_10 = [];
    dataForTransTarget = targetModel.refpts_10_5(idx_10_10,:);
    dataForTransMatch = matchModel.refpts_10_5(idx_10_10,:);
else
    disp('No scalp data type selected lol')
    return
end

[A,B] = DOTHUB_affineMap(dataForTransMatch,dataForTransTarget);
[i1,i2,i3] = ind2sub(size(matchModel.TissueMask),find(matchModel.TissueMask));
coords = [i1 i2 i3];
coords = DOTHUB_affineTrans(coords,A,B);
dummieImg = ones(size(matchModel.TissueMask));
[i1,i2,i3] = ind2sub(size(dummieImg),find(dummieImg));
dummieCoords = [i1 i2 i3];
[idx,~] = knnsearch(dummieCoords,coords);
newMask = zeros(size(dummieImg));
newMask(idx) = ones(length(idx),1);

newTissueMask = imfill(newMask,'holes');
dim = size(newTissueMask);
for i = 1:dim(2)-1
    slice = squeeze(newTissueMask(:,i,:));
    if max(slice(:)) == 0
        newTissueMask(:,i,:) = (newTissueMask(:,i+1,:));
    end
end
end