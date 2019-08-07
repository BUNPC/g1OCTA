%% g1-based OCT Angiogram, Repeat Ascan, GPU-based
% input: 
    % 1D array spectrum, nK*nXrpt*nX*nY, data format: ASCII int16
        % nK: spectrum pixel (camera elements); nXrpt: Ascan repeat;
        % nX: number of Ascans per Bscan; nY: number of Bscans for the whole volum
        % NOTE: the raw data for the whole volume is usually very large, it's recommended to process chunk by chunk
    % PRSinfo: processing information
    % PRSinfo.FWHM: Full width at Half Maxim, Amplitude, [transverse, axial], m
    % PRSinfo.Lam: [light source center, wavelength bandwidth], m
    % PRSinfo.Dim: [nz,nx,nyPchk,nTau]
    % PRSinfo.g1_Start: start time for g1 calculation
    % PRSinfo.g1_nt: number of time points for g1 calculation
    % PRSinfo.g1_ntau: number of g1 time lag
    % PRSinfo.intDk: OCT lambda to k space interpolation factor (calibration is required)
% subFunctions:
    % function [Dim, fNameBase, fIndex]=GetNameInfoRaw(filename0)
    % function DAT= ReadDat_int16(filePath, Dim, iseg, ARpt_extract,RptBscan) 
    % function RR = DAT2RR_GPU(Dat, intpDk)
    % function GG = RR2g1_GPU(RR, PRSinfo)
    % function g1AG=GG2g1AG_GPU(GG)
% output:
    % g1AG, [nz,nx,ny,2]
close all; clear;clc 
%% Set file location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datapath0  = 'H:\BU\PROJ - g1 OCTA\1205AnesIntralipidMouseDepth\ROI2\Dep1-Intralipid-DLSOCT';
[filename,datapath]=uigetfile(fullfile(datapath0,'*.dat'),'select file');
nameinfo=strsplit(filename,'-');
dx=1.5; % x pixel size, um
dz=3.3; % z pixel size, um
%% add MATLAB functions' path
addpath('.\SubFunctions') % sub functions
%% 0. data information
[Dim, fNameBase,fIndex]=GetNameInfoRaw(filename);
N_YChks=Dim.ny/10;   % default number of subYseg
%% 1. data process options
clear inputNum
prompt1={['Num Ysegs (nY=',num2str(Dim.ny),', nX=',num2str(Dim.nx),', NxRpt=',num2str(Dim.nxRpt),', NyRpt=',num2str(Dim.nyRpt),')'],...
    'RptA_Start (nARpt Process)','RptA_Inerval (nARpt Process)','RptA_n (nARpt Process)',...
    'g1_Start (g1 Process)','g1_nt (g1 Process)','g1_nTau (g1 Process)',...
    'N_Aline(DAQ)','N_Bscan(DAQ)','N_xRpt(DAQ)','N_yRpt(DAQ)',...
    'intDk','Aline rate (kHz)'};
inputNum1=inputdlg(prompt1,'', 1,{num2str(N_YChks),...
    '1','1',num2str(max(Dim.nxRpt,Dim.nyRpt)), ...
    '1','100','25', ...
    num2str(Dim.nx),num2str(Dim.ny),num2str(Dim.nxRpt),num2str(Dim.nyRpt),...
    '-0.43'});
N_YChks=str2num(inputNum1{1});    % number of chunks
% extract repeated alines from DAQ nxRpt, for repeated Ascan only
RptA_start=str2num(inputNum1{2});    % selecte start Aline repeat 
RptA_Interval=str2num(inputNum1{3});    % interval 
RptA_n=str2num(inputNum1{4});        % total number of extracted Alines
% g1 calculation parameters
PRSinfo.g1_Start=str2num(inputNum1{5});       % start time for g1 calculation
PRSinfo.g1_nt=str2num(inputNum1{6});       % number of time points for g1 calculation
PRSinfo.g1_ntau=str2num(inputNum1{7});       % number of g1 time lag
% DAQ info
Num_Aline=str2num(inputNum1{8});  % number of Aline
Num_Bscan=str2num(inputNum1{9});  % number of Bscans
n_xRpt=str2num(inputNum1{10});     % number of Ascan repeat
n_yRpt=str2num(inputNum1{11});     % number of Bscan repeat
PRSinfo.intDk=str2num(inputNum1{12});  % optional
%%%% 
ARpt_extract=[RptA_start,RptA_Interval,RptA_n];
%% 2. Select axial range for data processing %%%%%%
filePath=[datapath,filename];
disp(['Loading data to select the field of focus... ', num2str(floor(N_YChks/2)), ', ',datestr(now,'DD-HH:MM:SS')]);
DimChk=Dim;  DimChk.ny=1;
DAT = ReadDat_int16(filePath, DimChk, floor(N_YChks/2),ARpt_extract); % NII_ori: nk_Nx_ny,Nx=nt*nx;  floor(N_kfile/2)
RR=DAT2RR_GPU(DAT,PRSinfo.intDk);  % GPU
clear DAT
fig=figure;
imagesc(abs((squeeze(max(RR(:,:,:),[],3))))); caxis([0 5])
xlabel('Y');ylabel('Z');ylim([1 300]);title('MIP along X')
% caxis([0 5])
disp(['Select brain surface and stack start layer in figure']);
[XY_surf, Z_surf]=ginput(3);
close(fig);
prompt2={'Surface','Start Z_seg','End Z_seg'};
inputZseg=inputdlg(prompt2,'Z Segment parameter', 1,{num2str(floor(Z_surf(1))),num2str(floor(Z_surf(2))),num2str(floor(Z_surf(3)))});
z_seg_surf=str2num(inputZseg{1});
z_seg0=str2num(inputZseg{2});  % number of segments
LengthZ=str2num(inputZseg{3})-str2num(inputZseg{2});
%% 3. Load Data and process to g1, chunk by chunk %%%%%%%%%%%%%%%%%%
zRange=[z_seg0,z_seg0+LengthZ-1];
ARpt_extract=[RptA_start,RptA_Interval,RptA_n];
nyPerChk=floor(Dim.ny/N_YChks); % number of Bscans per ikfile
DimChk=Dim;
DimChk.ny=nyPerChk;
filePath=[datapath,filename];
GG=zeros(nz,nx,ny,PRSinfo.g1_ntau);
for iChk=1:N_YChks
    disp(['Start loading ith Chunk - ', num2str(iChk), ', ',datestr(now,'DD-HH:MM:SS')]);
    [DAT] = ReadDat_int16(filePath, DimChk, iChk, ARpt_extract); %, load the iChk chunk of data, nk_Nx_ny,Nx=nxRpt*nx 
    disp(['DAT2RR...',datestr(now,'DD-HH:MM:SS')]);
    RR0 = DAT2RR_GPU(DAT, PRSinfo.intDk);     % process raw spectrum data to spatial reflectivity data 
    RR=permute(reshape(RR0(zRange(1):zRange(2),:,:),[LengthZ,DimChk.nxRpt,DimChk.nx,DimChk.nyRpt,DimChk.ny]),[1 3 5 2 4]); % reshape RR from [nz,Nx,ny] to [nz,nx,ny,nxRpt,nyRpt], Nx=nx*nxRpt
    GG(:,:,(iChk-1)*nyPerChk+1:iChk*nyPerChk,:)=RR2g1_GPU(RR, PRSinfo); % calculate g1, GPU-based, ~40 times faster
end
%% 3.2 3x3x3 3D voxel average of g1
disp(['GG 3x3x3 3D averaging...',datestr(now,'DD-HH:MM:SS')]);
B_conv=ones(3,3,3);
for itau=1:PRSinfo.g1_ntau
    GG(:,:,:,itau)=convn(squeeze(GG(:,:,:,itau)),B_conv,'same')/numel(B_conv); % average neighboring 3X3X3, voxels
end
%% 3.3 GG to g1OCTA
disp(['GG to g1OCTA...',datestr(now,'DD-HH:MM:SS')]);
g1AG=GG2g1AG_GPU(GG); % GPU-based
save([datapath '\g1AG.mat'],'-v7.3','g1AG');
%% 4. g1OCTA result plot, MIP
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
g1AGMIP=squeeze(max(g1AG(:,:,:,1),[],1)).*squeeze(min(g1AG(:,:,:,2),[],1));
figure;
imagesc(g1AGMIP);colorbar; colormap (g1OCTAcmap);caxis([-1 1]); colorbar; axis equal; axis tight;axis off
saveas(gcf,[datapath,nameinfo{1},'-EnfaceMIP.jpg'],'jpg');
saveas(gcf,[datapath,nameinfo{1},'-EnfaceMIP.fig'],'fig');