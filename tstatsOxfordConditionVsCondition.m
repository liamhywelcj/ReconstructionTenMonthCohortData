function tstatImage = tstatsOxfordConditionVsCondition(twLower,twUpper,iChromo)
%This function produces the images of the contrast in responses between the two conditions

twLower = twLower + 20;
twUpper = twUpper + 20;

image_GM_Exl_TOTAL = [];
image_GM_Ctl_TOTAL = [];

for iSubject = 1:59
    dotimgFileName = ['~/OxfordReconstructionData/Participant',num2str(iSubject),'.dotimg'];
    if iChromo == 1
        load(dotimgFileName,'-mat','hbo')
        image_GM_Exl = hbo.gm(twLower:twUpper,:,3);
        image_GM_Ctl = hbo.gm(twLower:twUpper,:,2);
    else
        load(dotimgFileName,'-mat','hbr')
        image_GM_Exl = hbr.gm(twLower:twUpper,:,3);
        image_GM_Ctl = hbr.gm(twLower:twUpper,:,2);
    end
    image_GM_Exl_TOTAL = [image_GM_Exl_TOTAL; image_GM_Exl];
    image_GM_Ctl_TOTAL = [image_GM_Ctl_TOTAL; image_GM_Ctl];
    disp(iSubject)
end

load('OxfordArray.rmap','gmSurfaceMesh','-mat')
tVecContrast = zeros(length(gmSurfaceMesh.node),1);
pVecContrast = zeros(length(gmSurfaceMesh.node),1);
for i = 1:length(gmSurfaceMesh.node)
    [~,p,~,stats] = ttest(image_GM_Exl_TOTAL(:,i),image_GM_Ctl_TOTAL(:,i));
    pVecContrast(i) = p;
    tVecContrast(i) = stats.tstat;
end

% commonImageExl = mean(hbo_image_GM_Exl_TOTAL);
% commonImageCtl = mean(hbo_image_GM_Ctl_TOTAL);

a = 1;
for iSubject = 1:59
%     ID = subjectList(iSubject);
    load(['~/OxfordReconstructionData/Participant',num2str(iSubject),'.nirs'],'-mat','SD')
    if a == 1
        S = SD.MeasListAct;
    else
        S = S + SD.MeasListAct;
    end
    a = a+1;
    disp(a)
end
rmapFileName = '~OxfordArray.rmap';
load(rmapFileName,'-mat','headVolumeMesh')
basis = [40 40 40];
fineBasis = basis.*2;
disp('Mapping basis Jacobian to volume');
eltp = ones(length(headVolumeMesh.elem),1)*3;
hMesh = toastMesh(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4),eltp);
hBasis = toastBasis(hMesh,basis,fineBasis);
jacFileName = 'OxfordArray.jac';
load(jacFileName,'-mat','J')
J_GM_wav1 = J{1,1}.gm;
J_GM_wav2 = J{2,1}.gm;
J_basis_wav1 = J{1,1}.basis;
J_basis_wav2 = J{2,1}.basis;
disp('Mapping basis Jacobian to volume');
%Map back to volume
J_vol_wav1 = zeros(size(J_basis_wav1,1),length(headVolumeMesh.node));
J_vol_wav2 = zeros(size(J_basis_wav1,1),length(headVolumeMesh.node));
for i = 1:size(J_basis_wav1,1)
    J_vol_wav1(i,:) = hBasis.Map('S->M',J_basis_wav1(i,:));
    J_vol_wav2(i,:) = hBasis.Map('S->M',J_basis_wav2(i,:));
end

a = a-1;
thr = a*0.75;
S(S<thr) = 0;
S(S>=thr) = 1;

nirs = load('~/OxfordReconstructionData/Participant1.nirs','-mat');
SD = nirs.SD;
SD.MeasListAct = S;
[~, GMmask, ~, ~] = DOTHUB_MakeGMSensitivityMap(J_GM_wav1,J_GM_wav2,J_vol_wav1,J_vol_wav2,SD,0.01);
totalMask = GMmask;

pMaskContrast = pVecContrast<=0.01; 
pMaskContrast = double(pMaskContrast);
T_Contrast = pMaskContrast.*tVecContrast;
T_Contrast = T_Contrast.*totalMask';

tstatImage = T_Contrast;

end
