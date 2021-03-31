%% Head model creation
%Download the structural MRI data files from here: https://www.nitrc.org/projects/pediatricatlas
%This script relies on the iso2mesh toolbox, see http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download
MRISubj = 'infant-1yr-seg.nii';
seg = niftiread(MRISubj);
info = niftiinfo(MRISubj);
seg(seg == 250) = 3;
seg(seg == 150) = 2;
seg(seg == 10) = 1;
seg = double(seg);
headMask = niftiread('headMask.nii.gz');%preivously extracted head mask
headMask = double(headMask);
TissueMask = seg+headMask;

%% headVolumeMesh creation
maxVol = 1.25;voxelSize = 1;opt = 1;
dim = size(TissueMask);
[node,elem,face] = vol2mesh(uint8(mask),1:dim(1),1:dim(2),1:dim(3),opt,maxVol,1,'cgalmesh');
elem_or = meshreorient(node(:,1:3),elem(:,1:4)) ;
elem = [elem_or,elem(:,5)];
node(:,1:3) = node(:,1:3)*voxelSize;
headVolumeMesh.node = node;
headVolumeMesh.elem = elem;
headVolumeMesh.face = face;

mesh_orig = headVolumeMesh;
headVolumeMesh = DOTHUB_createNodalTissueInd(headVolumeMesh);
headVolumeMesh.elem = mesh_orig.elem;
headVolumeMesh.labels = {'ECT','CSF','GM','WM'};
[scalpSurf,~,~] = extractSurfacesFromMesh(headVolumeMesh);

%% gmSurfaceMesh creation
dim = size(TissueMask);
GM = TissueMask;
GM(GM<3) = 0;
GM = GM>0;
GM = double(GM);
GM_filled = zeros(dim);
for i = 1:dim(3)
    GM_filled(:,:,i) = imfill(GM(:,:,i),'holes');
end
dim = size(GM);
opt.radbound = 1;
[node_surf,face_surf] = vol2surf(GM_filled,1:dim(1),1:dim(2),1:dim(3),opt,1);
[node_r,face_r] = meshcheckrepair(node_surf,face_surf);
[newnode,newface] = meshresample(node_r,face_r,0.95);
[node_r,face_r] = meshcheckrepair(newnode,newface);
conn = meshconn(face_r,size(node_r,1));
node_s = smoothsurf(node_r,[],conn,10,0.95,'lowpass');
[node_s,face_r]=removeisolatednode(node_s,face_r);
[node_s,face_r] = meshcheckrepair(node_s,face_r);
[node_s,face_r] = surfreorient(node_s,face_r);
gmSurfaceMesh.node = node_s;
gmSurfaceMesh.face = face_r;
%% scalpSurfaceMesh creation
dim = size(headMask);
opt.radbound = 1;
opt.maxnode = 500000;
[node_surf,face_surf] = vol2surf(headMask,1:dim(1),1:dim(2),1:dim(3),opt,1);
[node_r,face_r] = meshcheckrepair(node_surf,face_surf);
[newnode,newface] = meshresample(node_r,face_r,0.95);
[node_r,face_r] = meshcheckrepair(newnode,newface);
conn = meshconn(face_r,size(node_r,1));
node_s = smoothsurf(node_r,[],conn,10,0.95,'lowpass');
[node_s,face_r]=removeisolatednode(node_s,face_r);
[node_s,face_r] = meshcheckrepair(node_s,face_r);
[node_s,face_r] = surfreorient(node_s,face_r);
scalpSurfaceMesh.node = node_s; 
scalpSurfaceMesh.face = face_r;

%% 10-5 positions computation
coordinates = [93 192 42;90 38 34;147 118 25;36 119 26;90 87 143];%reference points extracted from the MRI volume: Nz, Iz, Ar, Al and Cz
[landmarks,~,~] = DOTHUB_nearestNode(coordinates,scalpSurfaceMesh.node);
[refpts_10_5,~] = computeReferencePositions(scalpSurf(:,1:3),landmarks);
save('infantAtlasForReg.mshs','TissueMask','headVolumeMesh','gmSurfaceMesh','scalpSurfaceMesh','landmarks','refpts_10_5');
