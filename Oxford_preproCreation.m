function [prepro, preproFileName] = Oxford_preproCreation(nirsFileName,cfgFileName,nameSD3D) %,variableToReconstruct)

% ####################### OUTPUTS #########################################
%
% prepro            : The prepro stucture
%
% preproFileName    : The full pathname of the resulting prepro file
%
% ####################### Dependencies ####################################
% #########################################################################
% RJC, UCL, April 2020 and LCJ, UCL, February 2021

% MANAGE VARIABLES
% #########################################################################
[pathstr, name, ext] = fileparts(nirsFileName);
if isempty(ext) || ~strcmpi(ext,'.nirs')
    ext = '.nirs';
end
if isempty(pathstr) %Enforce full pathname
    pathstr = pwd;
end
nirsFileName = fullfile(pathstr,[name ext]);

% LOAD DATA
% #########################################################################
% Load .nirs
hmr = load(nirsFileName,'-mat');
load(nameSD3D,'-mat');hmr.SD = SD3D;hmr.SD3D = SD3D;

dcAvg = hmr.procResult.dcAvg;
dcAvgStd = hmr.procResult.dcAvgStd;
funcInd = find(strcmpi(hmr.procInput.procFunc.funcName,'hmrOD2Conc'));
DPFs = hmr.procInput.procFunc.funcParamVal{funcInd}{:};
dodRecon = DOTHUB_hmrHRFConc2OD(dcAvg,hmr.SD3D,DPFs);
tHRF = 1:1:length(dodRecon);
tRecon = tHRF;

% USE CODE SNIPPET FROM DOTHUB_writePREPRO to define filename and logData
[pathstr, name, ~] = fileparts(nirsFileName);
ds = datestr(now,'yyyymmDDHHMMSS');
preproFileName = fullfile(pathstr,[name '.prepro']);
logData(1,:) = {'Created on: '; ds};
logData(2,:) = {'Derived from data: ', nirsFileName};
logData(3,:) = {'Pre-processed using:', cfgFileName};

[prepro, preproFileName] = DOTHUB_writePREPRO(preproFileName,logData,dodRecon,tRecon,hmr.SD3D,hmr.s,dcAvg,dcAvgStd,tHRF);

end