function varargout = GUI_PLOT_g1OCTA(varargin)
% GUI_PLOT_g1OCTA MATLAB code for GUI_PLOT_g1OCTA.fig
%      GUI_PLOT_g1OCTA, by itself, creates a new GUI_PLOT_g1OCTA or raises the existing
%      singleton*.
%
%      H = GUI_PLOT_g1OCTA returns the handle to a new GUI_PLOT_g1OCTA or the handle to
%      the existing singleton*.
%
%      GUI_PLOT_g1OCTA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PLOT_g1OCTA.M with the given input arguments.
%
%      GUI_PLOT_g1OCTA('Property','Value',...) creates a new GUI_PLOT_g1OCTA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PLOT_g1OCTA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PLOT_g1OCTA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PLOT_g1OCTA

% Last Modified by GUIDE v2.5 17-Jul-2019 12:28:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PLOT_g1OCTA_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PLOT_g1OCTA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_PLOT_g1OCTA is made visible.
function GUI_PLOT_g1OCTA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PLOT_g1OCTA (see VARARGIN)
%% add MATLAB functions' path
handles.CodePath=pwd;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.defpath='H:';

handles.startZ=1;
handles.stackZ=100;
% Choose default command line output for GUI_PLOT_g1OCTA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_PLOT_g1OCTA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PLOT_g1OCTA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_LoadData.
function btn_LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to btn_LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
%% select file path %%%%%%%%%
defaultpath=handles.defpath;
[filename,datapath]=uigetfile(defaultpath);
handles.defpath=datapath;
handles.filename=filename;
%% load data %%%%%%%%%%
%%%%% input number of sub GG to be loaded %%%%%%%%
lding=msgbox(['Loading data...  ',datestr(now,'DD:HH:MM')]);
g1AG=LoadMAT(datapath,filename);
lded=msgbox(['Data loaded. ',datestr(now,'DD:HH:MM')]);
pause(1);
delete(lding); delete(lded);
[nz,nx,ny,nD]=size(g1AG);
%%%%%%%%%%%%%
prompt={'dX (um): ', 'dY(um): ', 'dZ size (um)'};
name='Enter Imaging info';
defaultvalue={'1.5','1.5','2.9'};
dXYZinput=inputdlg(prompt,name, 1, defaultvalue);
handles.Xcoor=[1:nx]*str2num(dXYZinput{1});
handles.Ycoor=[1:ny]*str2num(dXYZinput{2});
handles.Zcoor=[1:nz]*str2num(dXYZinput{3});

%% plot g1AG enface MIP
% handles.g1AGV=imgaussfilt3(g1AG(:,:,:,1),1.1);  % dynamic index  
% difAglGG=imgaussfilt3(g1AG(:,:,:,2),1.1);
% difAglGG(difAglGG>-3)=1;  % threshold to show only large descending flow
% handles.g1AGD=sign(difAglGG); % flow direction
handles.g1AGV=g1AG(:,:,:,1);
handles.g1AGD=g1AG(:,:,:,2);
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
vSign=squeeze(min(handles.g1AGD,[],1));
vSign(vSign~=-1)=1;
g1AGMIP=squeeze(max(handles.g1AGV,[],1)).*vSign;
axes(handles.axes1)
imagesc(g1AGMIP);
colorbar; colormap (g1OCTAcmap);caxis([-1 1]); 
axis equal; axis tight;
title('g1OCTA')

guidata(hObject, handles);


% --- Executes on button press in btn_plot.
function btn_plot_Callback(hObject, eventdata, handles)
% hObject    handle to btn_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[nz,nx,ny]=size(handles.g1AGV);

prompt={'MIP zStart', 'MIP zStack','Image refine','SideView(N:0;Y:1)'};
name='Enter Imaging info';
defaultvalue={'1',num2str(handles.startZ),num2str(min(nz,nz-handles.startZ)),'5'};
dXYZinput=inputdlg(prompt,name, 1, defaultvalue);
PlotSideView=str2num(dXYZinput{1});
zStart=str2num(dXYZinput{2});
zStack=str2num(dXYZinput{3});
rfn=str2num(dXYZinput{4});

handles.stackZ=zStack;
handles.startZ=zStart;
guidata(hObject, handles);  

%% PLOT
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
if PlotSideView==1 % plot SideView figures
    axes(handles.axes1);
    handles.slt(1)=str2num(get(handles.zStart,'string'));
    [handles.slt(3), handles.slt(2)]=ginput(1); % [y x]
    handles.slt=round(handles.slt);
    
    handles.fig=figure;
    set(handles.fig,'Position',[400 500 1400 800])
    subplot(2,2,1) % xz side view
    g1AGxy=squeeze(handles.g1AGV(handles.slt(1),:,:)).*squeeze(handles.g1AGD(handles.slt(1),:,:));
    g1AGxy(handles.slt(2)+[0 1],:)=1;
    g1AGxy(:,handles.slt(3)+[0 1])=1;
    imagesc(handles.Xcoor,handles.Ycoor,g1AGxy);
    colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
    axis equal; axis tight;
    title(['g1OCTA-XY, iz=',num2str(handles.slt(1))])
    xlabel('X [um]'); ylabel('Y [um]')
    
    subplot(2,2,2) % yz side view
    g1AGyz=squeeze(handles.g1AGV(:,handles.slt(2),:)).*squeeze(handles.g1AGD(:,handles.slt(2),:));
    g1AGyz(handles.slt(1),:)=1;
    imagesc(handles.Ycoor,handles.Zcoor,g1AGyz);
    colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
    axis equal; axis tight;
    title(['g1OCTA-XZ, iY=',num2str(handles.slt(2))])
    xlabel('X [um]'); ylabel('Z [um]')
    
    subplot(2,2,3) % XY single plan enface view
    g1AGMIP=squeeze(max(handles.g1AGV(zStart:zStart+zStack-1,:,:),[],1)).*squeeze(min(handles.g1AGD(zStart:zStart+zStack-1,:,:),[],1));
    g1AGMIP(handles.slt(2)+[0 1],:)=1;
    g1AGMIP(:,handles.slt(3)+[0 1])=1;
    imagesc(handles.Xcoor,handles.Ycoor,g1AGMIP);
    colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
    axis equal; axis tight;
    title('g1OCTA')
    xlabel('X [um]'); ylabel('Y [um]')
    
    subplot(2,2,4) % XY enface view MIP
    g1AGxz=squeeze(handles.g1AGV(:,:,handles.slt(3))).*squeeze(handles.g1AGD(:,:,handles.slt(3)));
    g1AGxz(handles.slt(1),:)=1;
    imagesc(handles.Xcoor,handles.Zcoor,g1AGxz);
    colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
    axis equal; axis tight;
    title(['g1OCTA-YZ, iX=',num2str(handles.slt(3))])
    xlabel('Y [um]'); ylabel('Z [um]')
else %% plot enface MIP only
    g1AGMIP=squeeze(max(handles.g1AGV(zStart:zStart+zStack-1,:,:),[],1)).*squeeze(min(handles.g1AGD(zStart:zStart+zStack-1,:,:),[],1));
    handles.fig=figure;
    imagesc(handles.Xcoor,handles.Ycoor,g1AGMIP);
    colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
    axis equal; axis tight;
    title('g1OCTA')
end
guidata(hObject, handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));

[nz,nx,ny] = size(handles.g1AGV);
set(hObject,'SliderStep',[1/(nz-1), 3/(nz-1)])
set(hObject,'Max',nz)
zStart=nz-min(round(get(hObject,'Value')),nz-1);
set(handles.zStart,'string',zStart);

g1AGMIP=squeeze(max(handles.g1AGV(zStart:min(zStart+zStack-1,nz),:,:),[],1)).*squeeze(min(handles.g1AGD(zStart:min(zStart+zStack-1,nz),:,:),[],1));
axes(handles.axes1)
imagesc(g1AGMIP);
colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
axis equal; axis tight;
title(['g1OCTA-XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
xlabel('X [pix]'); ylabel('Y [pix]')
guidata(hObject, handles);

function zStart_Callback(hObject, eventdata, handles)
% hObject    handle to zStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStart as text
%        str2double(get(hObject,'String')) returns contents of zStart as a double
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));
[nz,nx,ny]=size(handles.g1AGV);

g1AGMIP=squeeze(max(handles.g1AGV(zStart:min(zStart+zStack-1,nz),:,:),[],1)).*squeeze(min(handles.g1AGD(zStart:min(zStart+zStack-1,nz),:,:),[],1));
axes(handles.axes1)
imagesc(g1AGMIP);
colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
axis equal; axis tight;
title(['g1OCTA-XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
xlabel('X [pix]'); ylabel('Y [pix]')

function zStack_Callback(hObject, eventdata, handles)
% hObject    handle to zStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStack as text
%        str2double(get(hObject,'String')) returns contents of zStack as a double
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));
[nz,nx,ny]=size(handles.g1AGV);

g1AGMIP=squeeze(max(handles.g1AGV(zStart:min(zStart+zStack-1,nz),:,:),[],1)).*squeeze(min(handles.g1AGD(zStart:min(zStart+zStack-1,nz),:,:),[],1));
axes(handles.axes1)
imagesc(g1AGMIP);
colorbar; colormap (g1OCTAcmap);caxis([-1 1]);
axis equal; axis tight;
title(['g1OCTA-XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
xlabel('X [pix]'); ylabel('Y [pix]')

function btn_Save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
%% select file path %%%%%%%%%
defaultpath=handles.defpath;
disp('Saving data......')
saveas(handles.fig,[defaultpath,'AG-',handles.filename(1:end-4),'.fig']);
saveas(handles.fig,[defaultpath,'AG-',handles.filename(1:end-4),'.jpg']);
disp('Data saved!')



% --- Executes on button press in btn_reset.
function btn_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.defpath='H:';

handles.startZ=1;
handles.stackZ=100;
% Choose default command line output for GUI_PLOT_g1OCTA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function zStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function zStack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
