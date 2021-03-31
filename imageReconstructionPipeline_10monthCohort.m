%This script is the pipeline to reconstruct images from the fNIRS data from the 10-month cohort.
%This script relies on the DOT-HUB toolbox (https://github.com/liamhywelcj/DOT-HUB_toolbox) and the Homer2 toolbox (https://homer-fnirs.org/download/).
origMeshFileName = 'infantAtlasForReg.mshs';
load(origMeshFileName,'-mat','landmarks','refpts_10_5');
landmarksOrig = landmarks;
refpts_10_5_orig = refpts_10_5;
headCircumference = 45.4949; %cohort-mean head circumference
modelHeadCircumference = 4.737538345164647e+02; %previously calculated circumference of UNC 12-month model
scaleFactor = headCircumference/modelHeadCircumference; %the factor by which to scale the head model
reg.landmarks = landmarks*scaleFactor;
reg.headVolumeMesh = headVolumeMesh;
reg.headVolumeMesh.node(:,1:3) = reg.headVolumeMesh.node(:,1:3)*scaleFactor;
reg.scalpSurfaceMesh = scalpSurfaceMesh;
reg.scalpSurfaceMesh.node = reg.scalpSurfaceMesh.node*scaleFactor;
reg.refpts_10_5 = refpts_10_5*scaleFactor;
reg.gmSurfaceMesh = gmSurfaceMesh;
reg.gmSurfaceMesh.node = reg.gmSurfaceMesh.node*scaleFactor;

orig.TissueMask = TissueMask;
orig.landmarks = landmarks;
orig.refpts_10_5 = refpts_10_5;
reg.TissueMask = TissueMask;
reg.landmarks = landmarksOrig;
reg.refpts_10_5 = refpts_10_5_orig;
reg.headVolumeMesh = headVolumeMesh;
reg.scalpSurfaceMesh = scalpSurfaceMesh;
reg.gmSurfaceMesh = gmSurfaceMesh;

TissueMaskNew = affineRegisterTissueMask(orig,reg,'landmarks');
[A,B] = DOTHUB_affineMap(orig.landmarks,reg.landmarks);
reg.gmSurfaceMesh.node = DOTHUB_affineTrans(reg.gmSurfaceMesh.node,A,B);
reg.scalpSurfaceMesh.node = DOTHUB_affineTrans(reg.scalpSurfaceMesh.node,A,B);

save('infantAtlasWarped.mshs','-struct','reg')
SpringRelaxationRegistration; %this is a function
nameSD3D = 'OxfordArray.SD3D';
[rmap, rmapFileName] = OxfordSaveRMAP(nameSD3D,origMeshFileName);
[jac, jacFileName] = DOTHUB_makeToastJacobian(rmapFileName,[40 40 40]);

for iSubject = 1:59
    nirsFileName = ['Participant',num2str(iSubject),'.nirs'];
    [prepro, preproFileName] = Oxford_preproCreation(nirsFileName,'For reconstruction',nameSD3D);
    [dotimg,~] = DOTHUB_reconstruction(preproFileName,jacFileName,[],rmapFileName,'reconSpace','cortex','saveVolumeImages',false,'imagetype','both');
end