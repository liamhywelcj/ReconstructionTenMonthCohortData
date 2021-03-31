function pos = SpringRelaxationRegistration
%this relies on Homer2 functions, see https://homer-fnirs.org/download/
clear variables

load('arrayWithSprings.SD','-mat','SD'); 
load('infantAtlasWarped.mshs','-mat','refpts_10_5','TissueMask')
refpts.pos = refpts_10_5;
refpts.labels = {'Cz';'Fpz';'AFpz';'AFz';'AFFz';'Fz';'FFCz';'FCz';'FCCz';'CCPz';'CPz';'CPPz';'Pz';'PPOz';'POz';'POOz';'Oz';...
    'T7';'T7h';'C5';'C5h';'C3';'C3h';'C1';'C1h';'C2h';'C2';'C4h';'C4';'C6h';'C6';'T8h';'T8';'Fp1h';'Fp1';'AFp7';...
    'AF7';'AFF7';'F7';'FFT7';'FT7';'FTT7';'TTP7';'TP7';'TPP7';'P7';'PPO7';'PO7';'POO7';'O1';'O1h';'Fp2h';'Fp2';...
    'AFp8';'AF8';'AFF8';'F8';'FFT8';'FT8';'FTT8';'TTP8';'TP8';'TPP8';'P8';'PPO8';'PO8';'POO8';'O2';'O2h';'F7h';...
    'F5';'F5h';'F3';'F3h';'F1';'F1h';'F2h';'F2';'F4h';'F4';'F6h';'F6';'F8h';'P7h';'P5';'P5h';'P3';'P3h';'P1';...
    'P1h';'P2h';'P2';'P4h';'P4';'P6h';'P6';'P8h';'AF7h';'AF5';'AF5h';'AF3';'AF3h';'AF1';'AF1h';'AF2h';'AF2';...
    'AF4h';'AF4';'AF6h';'AF6';'AF8h';'PO7h';'PO5';'PO5h';'PO3';'PO3h';'PO1';'PO1h';'PO2h';'PO2';'PO4h';'PO4';...
    'PO6h';'PO6';'PO8h';'FT7h';'FC5';'FC5h';'FC3';'FC3h';'FC1';'FC1h';'FC2h';'FC2';'FC4h';'FC4';'FC6h';'FC6';...
    'FT8h';'TP7h';'CP5';'CP5h';'CP3';'CP3h';'CP1';'CP1h';'CP2h';'CP2';'CP4h';'CP4';'CP6h';'CP6';'TP8h';'FTT7h';...
    'FCC5';'FCC5h';'FCC3';'FCC3h';'FCC1';'FCC1h';'FCC2h';'FCC2';'FCC4h';'FCC4';'FCC6h';'FCC6';'FTT8h';'TTP7h';'CCP5';...
    'CCP5h';'CCP3';'CCP3h';'CCP1';'CCP1h';'CCP2h';'CCP2';'CCP4h';'CCP4';'CCP6h';'CCP6';'TTP8h';'FFT7h';'FFC5';'FFC5h';...
    'FFC3';'FFC3h';'FFC1';'FFC1h';'FFC2h';'FFC2';'FFC4h';'FFC4';'FFC6h';'FFC6';'FFT8h';'TPP7h';'CPP5';'CPP5h';'CPP3';'CPP3h';...
    'CPP1';'CPP1h';'CPP2h';'CPP2';'CPP4h';'CPP4';'CPP6h';'CPP6';'TPP8h';'AFF7h';'AFF5';'AFF5h';'AFF3';'AFF3h';'AFF1';'AFF1h';'AFF2h';...
    'AFF2';'AFF4h';'AFF4';'AFF6h';'AFF6';'AFF8h';'PPO7h';'PPO5';'PPO5h';'PPO3';'PPO3h';'PPO1';'PPO1h';'PPO2h';'PPO2';'PPO4h';'PPO4';...
    'PPO6h';'PPO6';'PPO8h';'AFp5';'AFp3';'AFp1';'AFp2';'AFp4';'AFp6';'POO5';'POO3';'POO1';'POO2';'POO4';'POO6'};

probe.optpos = [SD.SrcPos; SD.DetPos; SD.DummyPos];
probe.sl = SD.SpringList;
probe.al = SD.AnchorList;
mask = TissueMask; 
mask(mask>0) = 1;
imgData = mask;
dims = size(imgData);
headvol.img = imgData;
headvol.center = round(dims/2);

[optconn, ~] = spring2posprobe(probe, refpts, headvol); 
anchor_points(:,1) = [16,33:42,18,21,27,32];
anchor_points(:,2:4) = refpts.pos([2;1;6;22;29;24;27;39;57;18;33;89;92;35;53],:);% polhemus_anchor_points = refpts.pos([1;13;6],:);
posprobe_data = gen_positionprobe_dat(probe.optpos, optconn, anchor_points);  
if isempty(posprobe_data)
    return;
end

permuteAxes = [1     2     3];
scale = 1;
poso = [0 0 0];
posAttract = headvol.center;
nIter = 1000;
[pos,~] = positionprobe(posprobe_data,permuteAxes,scale,poso,imgData,posAttract,nIter);
idxSrc = 1:16;
idxDet = 17:32;
srcPos = pos(idxSrc,:);
detPos = pos(idxDet,:);

[srcPos,~,~] = DOTHUB_nearestNode(srcPos,scalpSurfaceMesh.node);
[detPos,~,~] = DOTHUB_nearestNode(detPos,scalpSurfaceMesh.node);

P = anchor_points(:,2:4);
P = bring_pts_to_surf_SB(scalpSurfaceMesh.node,P);

figure;plotmesh(scalpSurfaceMesh.node,scalpSurfaceMesh.face)
hold on
plotmesh(P,'g.','markersize',30)
plotmesh(srcPos,'r.','markersize',30)
plotmesh(detPos,'b.','markersize',30)
set(gcf,'color','w')

SD3D = SD;
SD3D.SrcPos = srcPos;
SD3D.DetPos = detPos;
SD3D.nSrcs = length(idxSrc);
SD3D.nDets = length(idxDet);
SD3D.Landmarks = landmarks;
SD3D.MeasListAct = ones(length(SD.MeasList),1);
nameSD3D = 'OxfordArray.SD3D';
save(nameSD3D,'SD3D')

end



