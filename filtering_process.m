function varargout = filtering_process(varargin)
%FILTERING_PROCESS MATLAB code file for filtering_process.fig
%      FILTERING_PROCESS, by itself, creates a new FILTERING_PROCESS or raises the existing
%      singleton*.
%
%      H = FILTERING_PROCESS returns the handle to a new FILTERING_PROCESS or the handle to
%      the existing singleton*.
%
%      FILTERING_PROCESS('Property','Value',...) creates a new FILTERING_PROCESS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to filtering_process_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      FILTERING_PROCESS('CALLBACK') and FILTERING_PROCESS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in FILTERING_PROCESS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help filtering_process

% Last Modified by GUIDE v2.5 08-Jun-2022 19:56:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @filtering_process_OpeningFcn, ...
                   'gui_OutputFcn',  @filtering_process_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before filtering_process is made visible.
function filtering_process_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for filtering_process
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% fs=varargin{1};
% x=varargin{2};
% t=varargin{3};
% UIWAIT makes filtering_process wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = filtering_process_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function f_pass_Callback(hObject, eventdata, handles)
% hObject    handle to f_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_pass as text
%        str2double(get(hObject,'String')) returns contents of f_pass as a double


% --- Executes during object creation, after setting all properties.
function f_pass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_cut_Callback(hObject, eventdata, handles)
% hObject    handle to f_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_cut as text
%        str2double(get(hObject,'String')) returns contents of f_cut as a double


% --- Executes during object creation, after setting all properties.
function f_cut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function decay_max_Callback(hObject, eventdata, handles)
% hObject    handle to decay_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of decay_max as text
%        str2double(get(hObject,'String')) returns contents of decay_max as a double


% --- Executes during object creation, after setting all properties.
function decay_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to decay_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function decay_min_Callback(hObject, eventdata, handles)
% hObject    handle to decay_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of decay_min as text
%        str2double(get(hObject,'String')) returns contents of decay_min as a double


% --- Executes during object creation, after setting all properties.
function decay_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to decay_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filtering_process.
function filtering_process_Callback(hObject, eventdata, handles)
% hObject    handle to filtering_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles_part1=guihandles(part1);
[fs,x,t,s]=part1;
IIR=get(handles.IIR,'Value');
FIR=get(handles.FIR,'Value');
f_pass1= str2double(get(handles.f_pass,'String'));% 通带频率
f_cut1= str2double(get(handles.f_cut,'String'));%截至频率f_pass= str2double(get(handles.f_pass,'String'));% 通带频率
f_pass2= str2double(get(handles.f_pass2,'String'));% 通带频率
f_cut2= str2double(get(handles.f_cut2,'String'));%截至频率
decay_max= str2double(get(handles.decay_max,'String'));
decay_min= str2double(get(handles.decay_min,'String'));
wp = [f_pass1 f_pass2]*2/fs;%将模拟频率变换到数字域w=2pi*fp/fs [0.1pi-0.2pi]不乘Pi,单位为pi，通带截止频率
ws = [f_cut1 f_cut2]*2/fs;%阻带起始频率
rp=decay_max;
rs=decay_min;
plot(handles.axes1,t,x(1,:));
xlabel(handles.axes1,'接收信号的时间/s');
title(handles.axes1,'阵元1接收的信号波形');

if IIR==1

[N,wn]=buttord(wp,ws,rp,rs);%N为最小阶数，Wn为3db带宽
[b1,a1]=butter(N,wn);
%w=linspace(0,pi,200);
[h,f]=freqz(b1,a1,8000,fs);
plot(handles.axes4,f,10*log(abs(h)));
ylim(handles.axes4,[-50,0])
xlabel(handles.axes4,'频率/Hz');
ylabel(handles.axes4,'幅度/dB');
title(handles.axes4,'IIR巴特沃斯滤波器');

%plot(handles.axes1,h)
y_iir=filter(b1,a1,x(1,:));%滤波
plot(handles.axes2,t,y_iir);
ylabel(handles.axes2,'x')
xlabel(handles.axes2,'接收信号的时间/s')
title(handles.axes2,'阵元1接收的信号波形（滤波后）');

Y_iir_fft=fft(y_iir);%FFT
L1=length(Y_iir_fft);
f1=(0:L1/2-1)*fs/L1;
plot(handles.axes3,f1,abs(Y_iir_fft(1,1:L1/2)));
ylabel(handles.axes3,'Y')
xlabel(handles.axes3,'频率/Hz')
title(handles.axes3,'阵元1接收的信号频域（滤波后）');
else
    [N,wn]=buttord(wp,ws,rp,rs);
    b=fir1(N,wn);%默认海明窗 N为滤波器阶数，越大主瓣越窄，过渡带越窄N=64
    a=1;
    y_fir=filter(b,a,x(1,:));%a=1幅值可以到0db
    [h,f]=freqz(b,a,8000,fs);
plot(handles.axes4,f,10*log(abs(h)));
ylim(handles.axes4,[-50,0])
xlabel(handles.axes4,'频率/Hz');
ylabel(handles.axes4,'幅度/dB');
    title(handles.axes4,'FIR滤波器');

plot(handles.axes2,t,y_fir);
ylabel(handles.axes2,'x')
xlabel(handles.axes2,'接收信号的时间/s')
title(handles.axes2,'阵元1接收的信号波形（滤波后）');

Y_fir_fft=fft(y_fir);%FFT
L1=length(Y_fir_fft);
f1=(0:L1/2-1)*fs/L1;

plot(handles.axes3,f1,abs(Y_fir_fft(1,1:L1/2)));
ylabel(handles.axes3,'Y')
xlabel(handles.axes3,'频率/Hz')
title(handles.axes3,'阵元1接收的信号频域（滤波后）');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); %指定需要清空的坐标轴
cla reset;
axes(handles.axes2); %指定需要清空的坐标轴
cla reset;
axes(handles.axes3); %指定需要清空的坐标轴
cla reset;
axes(handles.axes4); %指定需要清空的坐标轴
cla reset;

handles.waveform.Toolbar.Visible = 'on';



function f_pass2_Callback(hObject, eventdata, handles)
% hObject    handle to f_pass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_pass2 as text
%        str2double(get(hObject,'String')) returns contents of f_pass2 as a double


% --- Executes during object creation, after setting all properties.
function f_pass2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_pass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_cut2_Callback(hObject, eventdata, handles)
% hObject    handle to f_cut2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_cut2 as text
%        str2double(get(hObject,'String')) returns contents of f_cut2 as a double


% --- Executes during object creation, after setting all properties.
function f_cut2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_cut2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save1_Callback(hObject, eventdata, handles)
% hObject    handle to save1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFig = figure;%由于直接保存axes1上的图像有困难，所以保存在新建的figure中的谱图
set(newFig,'Visible','off')%设置新建的figure为不可见
newAxes = copyobj(handles.axes1,newFig);   %将axes1中的图复制到新建的figure中
set(newAxes,'Units','default','Position','default');    % 设置图显示的位置
[filename,pathname] = uiputfile({ '*.jpg','figure type(*.jpg)'}, '保存原始波形');
if isequal(filename,0)||isequal(pathname,0)%如果用户选择“取消”，则退出
    return;
else
    fpath=fullfile(pathname,filename);
end
print(newFig,'-djpeg','-r600',fpath);

% --------------------------------------------------------------------
function save2_Callback(hObject, eventdata, handles)
% hObject    handle to save2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFig = figure;%由于直接保存axes1上的图像有困难，所以保存在新建的figure中的谱图
set(newFig,'Visible','off')%设置新建的figure为不可见
newAxes = copyobj(handles.axes2,newFig);   %将axes1中的图复制到新建的figure中
set(newAxes,'Units','default','Position','default');    % 设置图显示的位置
[filename,pathname] = uiputfile({ '*.jpg','figure type(*.jpg)'}, '保存原始波形');
if isequal(filename,0)||isequal(pathname,0)%如果用户选择“取消”，则退出
    return;
else
    fpath=fullfile(pathname,filename);
end
print(newFig,'-djpeg','-r600',fpath);

% --------------------------------------------------------------------
function save3_Callback(hObject, eventdata, handles)
% hObject    handle to save3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFig = figure;%由于直接保存axes1上的图像有困难，所以保存在新建的figure中的谱图
set(newFig,'Visible','off')%设置新建的figure为不可见
newAxes = copyobj(handles.axes3,newFig);   %将axes1中的图复制到新建的figure中
set(newAxes,'Units','default','Position','default');    % 设置图显示的位置
[filename,pathname] = uiputfile({ '*.jpg','figure type(*.jpg)'}, '保存原始波形');
if isequal(filename,0)||isequal(pathname,0)%如果用户选择“取消”，则退出
    return;
else
    fpath=fullfile(pathname,filename);
end
print(newFig,'-djpeg','-r600',fpath);

% --------------------------------------------------------------------
function save4_Callback(hObject, eventdata, handles)
% hObject    handle to save4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newFig = figure;%由于直接保存axes1上的图像有困难，所以保存在新建的figure中的谱图
set(newFig,'Visible','off')%设置新建的figure为不可见
newAxes = copyobj(handles.axes4,newFig);   %将axes1中的图复制到新建的figure中
set(newAxes,'Units','default','Position','default');    % 设置图显示的位置
[filename,pathname] = uiputfile({ '*.jpg','figure type(*.jpg)'}, '保存原始波形');
if isequal(filename,0)||isequal(pathname,0)%如果用户选择“取消”，则退出
    return;
else
    fpath=fullfile(pathname,filename);
end
print(newFig,'-djpeg','-r600',fpath);
