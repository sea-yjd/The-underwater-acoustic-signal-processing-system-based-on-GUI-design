function varargout = part1(varargin)
% PART1 MATLAB code for part1.fig
%      PART1, by itself, creates a new PART1 or raises the existing
%      singleton*.
%
%      H = PART1 returns the handle to a new PART1 or the handle to
%      the existing singleton*.
%
%      PART1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PART1.M with the given input arguments.
%
%      PART1('Property','Value',...) creates a new PART1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before part1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to part1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help part1

% Last Modified by GUIDE v2.5 08-Jun-2022 19:14:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @part1_OpeningFcn, ...
                   'gui_OutputFcn',  @part1_OutputFcn, ...
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


% --- Executes just before part1 is made visible.
function part1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to part1 (see VARARGIN)

% Choose default command line output for part1
global fs;
global t;
global x;
global s;
global f0;
global Bw;
global Pw;
global t0;
global theata0;
global num;
global s_s;
global LFM;
global CW;
global SNR;
handles.output =fs;
handles.second_output = x;
handles.third_output = t;
handles.forth_output = s;
handles.fifth_output = f0;
handles.sixth_output = Bw;
handles.seventh_output = Pw;
handles.eigth_output = t0;
handles.ninth_output = theata0;
handles.tenth_output = num;
handles.eleventh_output = s_s;
handles.twelfth_output = LFM;
handles.thirteenth_output = CW;
handles.forteenth_output = SNR;
% handles.second_output = x;
% handles.third_output = t;
% Update handles structure
guidata(hObject, handles);

%filteringHandles = varargin{1};

% UIWAIT makes part1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = part1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.second_output;
varargout{3} = handles.third_output;
varargout{4} = handles.forth_output;
varargout{5} = handles.fifth_output;
varargout{6} = handles.sixth_output;
varargout{7} = handles.seventh_output;
varargout{8} = handles.eigth_output;
varargout{9} = handles.ninth_output;
varargout{10} = handles.tenth_output;
varargout{11} = handles.eleventh_output;
varargout{12} = handles.twelfth_output;
varargout{13} = handles.thirteenth_output;
varargout{14} = handles.forteenth_output;



% --- Executes on button press in sig_gen.
function sig_gen_Callback(hObject, eventdata, handles)
% hObject    handle to sig_gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fs;
global t;
global x;
global s;
global f0;
global Bw;
global Pw;
global t0;
global theata0;
global num;
global s_s;
global LFM;
global CW;
global SNR;
f0=get(handles.f0,'String');f0=str2double(f0);
LFM=get(handles.LFM,'Value');%LFM=str2double(LFM);
CW=get(handles.CW,'Value');%CW=str2double(CW);
Bw=get(handles.Bw,'String');Bw=str2double(Bw);
fs=get(handles.fs,'String');fs=str2double(fs);
Pw=get(handles.Pw,'String');Pw=str2double(Pw);
Pt=get(handles.Pt,'String');Pt=str2double(Pt);
SNR=get(handles.SNR,'String');SNR=str2double(SNR);
t0=get(handles.t0,'String');t0=str2double(t0);
theata0=get(handles.theata0,'String');theata0=str2double(theata0);

M=Bw/Pw;
num=get(handles.num,'String');num=str2double(num);
d=get(handles.d,'String');d=str2double(d);
l=length(num2str(fs));
delta=d*cosd(theata0);
c=1500;
tao=round(delta/c,l-1);

hold on;

for k=1:num
    t1=0:1/fs:t0+(k-1)*tao-1/fs;
    t2=t0+(k-1)*tao:1/fs:t0+Pw+(k-1)*tao-1/fs;
    t3=t0+Pw+(k-1)*tao:1/fs:t0+Pt-Pw-1/fs;
	t=[t1,t2,t3];
    s1=zeros(1,length(t1));
    sig_LFM=cos(2*pi*((f0-Bw/2-M)*t2+0.5*M*t2.^2));
    sig_CW=1/sqrt(Pw)*cos(2*pi*f0*t2);
    t_1=0:1/fs:Pw-1/fs;
    if LFM==1
         s2=sig_LFM;
         s_s=cos(2*pi*((f0-Bw/2-M)*t_1+0.5*M*t_1.^2));
    else
         s2=sig_CW;
         s_s=1/sqrt(Pw)*cos(2*pi*f0*t_1);
    end
    s3=zeros(1,length(t3));
    s(k,:)=[s1,s2,s3];
    x(k,:)=awgn(s(k,:),SNR,'measured');
    plot(handles.multi_wave,t,x(k,:));
    xlabel(handles.multi_wave,'接收信号的时间/s');
    ylabel(handles.multi_wave,'信号幅值');
    legend(handles.multi_wave,'location','NorthWest');
    title(handles.multi_wave,['多个阵元接收信号的时域波形']);
end
hold off;
plot(handles.waveform,t,x(1,:));
xlabel(handles.waveform,'接收信号的时间/s');
title(handles.waveform,'阵元1接收的主动信号波形');



% --- Executes on button press in sig_fft.
function sig_fft_Callback(hObject, eventdata, handles)
% hObject    handle to sig_fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t;
global x;
global fs;
global s;
% f0=get(handles.f0,'String');f0=str2double(f0);
% LFM=get(handles.LFM,'Value');%LFM=str2double(LFM);
% CW=get(handles.CW,'Value');%CW=str2double(CW);
% Bw=get(handles.Bw,'String');Bw=str2double(Bw);
% fs=get(handles.fs,'String');fs=str2double(fs);
% Pw=get(handles.Pw,'String');Pw=str2double(Pw);
% Pt=get(handles.Pt,'String');Pt=str2double(Pt);
% SNR=get(handles.SNR,'String');SNR=str2double(SNR);
% t0=get(handles.t0,'String');t0=str2double(t0);
% theata0=get(handles.theata0,'String');theata0=str2double(theata0);
% 
% M=Bw/Pw;
% num=get(handles.num,'String');num=str2double(num);
% d=get(handles.d,'String');d=str2double(d);
% l=length(num2str(fs));
% delta=d*cosd(theata0);
% c=1500;
% tao=round(delta/c,l-1);
%     k=1;
%     t1=0:1/fs:t0+(k-1)*tao-1/fs;
%     t2=t0+(k-1)*tao:1/fs:t0+Pw+(k-1)*tao-1/fs;
%     t3=t0+Pw+(k-1)*tao:1/fs:t0+Pt-Pw-1/fs;
% 	t=[t1,t2,t3];
%     s1=zeros(1,length(t1));
%     sig_LFM=cos(2*pi*((f0-Bw/2-M)*t2+0.5*M*t2.^2));
%     sig_CW=1/sqrt(Pw)*cos(2*pi*f0*t2);
%     if LFM==1
%          s2=sig_LFM;
%          
%     else
%          s2=sig_CW;
%     end
%     s3=zeros(1,length(t3));
%     s(k,:)=[s1,s2,s3];
%     x(k,:)=awgn(s(k,:),SNR,'measured');
sig_fft=fft(x(1,:));
L1=length(sig_fft);
f1=(0:L1/2-1)*fs/L1;
plot(handles.f_wave,f1,abs(sig_fft(1,1:L1/2)));
xlabel(handles.f_wave,'频率/Hz');
title(handles.f_wave,'阵元1主动信号频域(从接收信号开始）');



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f0_Callback(hObject, eventdata, handles)
% hObject    handle to f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0 as text
%        str2double(get(hObject,'String')) returns contents of f0 as a double


% --- Executes during object creation, after setting all properties.
function f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bw_Callback(hObject, eventdata, handles)
% hObject    handle to Bw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bw as text
%        str2double(get(hObject,'String')) returns contents of Bw as a double


% --- Executes during object creation, after setting all properties.
function Bw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_Callback(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs as text
%        str2double(get(hObject,'String')) returns contents of fs as a double


% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pw_Callback(hObject, eventdata, handles)
% hObject    handle to Pw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pw as text
%        str2double(get(hObject,'String')) returns contents of Pw as a double


% --- Executes during object creation, after setting all properties.
function Pw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pt_Callback(hObject, eventdata, handles)
% hObject    handle to Pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pt as text
%        str2double(get(hObject,'String')) returns contents of Pt as a double


% --- Executes during object creation, after setting all properties.
function Pt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SNR_Callback(hObject, eventdata, handles)
% hObject    handle to SNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR as text
%        str2double(get(hObject,'String')) returns contents of SNR as a double


% --- Executes during object creation, after setting all properties.
function SNR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t0_Callback(hObject, eventdata, handles)
% hObject    handle to t0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t0 as text
%        str2double(get(hObject,'String')) returns contents of t0 as a double


% --- Executes during object creation, after setting all properties.
function t0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theata0_Callback(hObject, eventdata, handles)
% hObject    handle to theata0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theata0 as text
%        str2double(get(hObject,'String')) returns contents of theata0 as a double


% --- Executes during object creation, after setting all properties.
function theata0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theata0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CW.
function CW_Callback(hObject, eventdata, handles)
% hObject    handle to CW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CW


% --- Executes on button press in LFM.
function LFM_Callback(hObject, eventdata, handles)
% hObject    handle to LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LFM



function num_Callback(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num as text
%        str2double(get(hObject,'String')) returns contents of num as a double


% --- Executes during object creation, after setting all properties.
function num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_Callback(hObject, eventdata, handles)
% hObject    handle to d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d as text
%        str2double(get(hObject,'String')) returns contents of d as a double


% --- Executes during object creation, after setting all properties.
function d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in waveform_selection.
function waveform_selection_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in waveform_selection 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.waveform); %指定需要清空的坐标轴
cla reset;
axes(handles.f_wave); %指定需要清空的坐标轴
cla reset;
axes(handles.multi_wave); %指定需要清空的坐标轴
cla reset;

handles.waveform.Toolbar.Visible = 'on';


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x;
global t;
global fs;
global s;
global f0;
global num;
global theata0;
global SNR;
c = 1500;
num = str2double(get(handles.num,'String'));% 阵元数目
d = str2double(get(handles.d,'String'));  % 阵列间距
f0 = str2double(get(handles.line_f,'String')); % 线谱频率
B = str2double(get(handles.band_width,'String'));  % 宽带谱带宽
SNR = str2double(get(handles.SNR_acc,'String'));
fs=str2double(get(handles.fs,'String'));
Pw=get(handles.Pw,'String');Pw=str2double(Pw);
Pt=get(handles.Pt,'String');Pt=str2double(Pt);
t0=get(handles.t0,'String');t0=str2double(t0);
theata0 = str2double(get(handles.theta_object,'String'));

T = t0+Pt-Pw-1/fs;  % 接收信号时长
t = 0: 1/fs: T;
% 阵元上的其他噪声
noise_env = randn(num,length(t));

s1_amplitude = randn(1);  % 线谱幅度
s1_phase = rand(1) * 2 * pi;  % 线谱相位
s1_t = s1_amplitude * cos(2*pi*f0*t+s1_phase);  % 线谱成分
% 宽带噪声
fl = f0 - B/2;
fh = f0 + B/2;
s2_t = randn(1,length(t));
s2_fft = fft(s2_t);
num_fl = floor(fl/fs*length(t));
num_fh = ceil(fh/fs*length(t));
s2_fft = [zeros(1,num_fl-1) s2_fft(num_fl:num_fh) zeros(1,length(t)-num_fh)];  
s2_t = real(ifft(s2_fft)); 
%海洋环境噪声
s_t = s1_t + s2_t;
power_s = 10^(SNR/10) * var(noise_env(1,:));              
s = sqrt(power_s) / std(s_t) * s_t;   
x = s + noise_env(1,:);
% 第一个阵元信号时域波形
plot(handles.waveform,t,x);
title(handles.waveform,'第一个阵元信号时域波形');
xlabel(handles.waveform,'接收信号的时间/s');
xlim(handles.waveform,[0,T]);
ylim(handles.waveform,[min(x),max(x)])

% 第一个阵元信号频谱
x_fft = fft(x);
L1=length(x_fft);
f1=(0:length(x_fft)/2-1)/length(x_fft)*fs;
plot(handles.f_wave,f1,abs(x_fft(1,1:L1/2)),'LineWidth',1.0)
title(handles.f_wave,'第一个阵元信号频谱')
xlabel(handles.f_wave,'频率/Hz')
xlim(handles.f_wave,[0,fs/2])

%% 阵元上的信号
delaynum_object = round((0:num-1).'*d*cos(theata0)/c*fs);  % 以第一个阵元上的信号为参考信号   %% 
v_object = exp(-1j*2*pi/length(t)*delaynum_object*(0:length(t)-1));  
s_fft = fft(s_t);
s_fft_vector = v_object .* s_fft;
s_t_vector = ifft(s_fft_vector.').';
x_t_vector = s_t_vector + randn(num,length(t));
x_fft_vector = fft(x_t_vector.').';
%% ULA输出信号%感觉不太对 噪声衰减的太多了 线谱也被衰减了
theta_ULA = 90; 
delaynum_array = round((0:num-1).'*d*cos(theta_ULA)/c*fs);
weight_amplitude = ones(num,1) / num ;%均匀加权
weight_phase = exp(-1j*2*pi/length(t)*delaynum_array*(0:length(t)-1));
weight = weight_amplitude .* weight_phase;
y_fft = sum(conj(weight) .* x_fft_vector,1);
y_t = ifft(y_fft);
% ULA输出信号的时域波形

plot(handles.multi_wave,t,real(y_t),'linewidth',1.0);
title(handles.multi_wave,'ULA均匀加权后输出信号信号时域波形');
xlim(handles.multi_wave,[0,T]);
ylim(handles.multi_wave,[min(x),max(x)])
xlabel(handles.multi_wave,'时间/s');

% ULA输出信号的频谱
% 
% plot(handles.ULA,(1:length(y_fft))/length(y_fft)*fs,abs(y_fft)/max(abs(y_fft)),'LineWidth',1.0)
% title(handles.ULA,'ULA均匀加权后输出信号信号频谱');
% xlabel(handles.ULA,'频率(Hz)')
% xlim(handles.ULA,[0,fs/2])



function line_f_Callback(hObject, eventdata, handles)
% hObject    handle to line_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of line_f as text
%        str2double(get(hObject,'String')) returns contents of line_f as a double


% --- Executes during object creation, after setting all properties.
function line_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to line_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function band_width_Callback(hObject, eventdata, handles)
% hObject    handle to band_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of band_width as text
%        str2double(get(hObject,'String')) returns contents of band_width as a double


% --- Executes during object creation, after setting all properties.
function band_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to band_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SNR_acc_Callback(hObject, eventdata, handles)
% hObject    handle to SNR_acc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SNR_acc as text
%        str2double(get(hObject,'String')) returns contents of SNR_acc as a double


% --- Executes during object creation, after setting all properties.
function SNR_acc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SNR_acc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theta_object_Callback(hObject, eventdata, handles)
% hObject    handle to theta_object (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta_object as text
%        str2double(get(hObject,'String')) returns contents of theta_object as a double


% --- Executes during object creation, after setting all properties.
function theta_object_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_object (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g=filtering_process;
set(g,'Visible','on');
fig=openfig('filtering_process.fig');
handles = guihandles(fig);
guidata(fig, handles);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in spect.
function spect_Callback(hObject, eventdata, handles)
% hObject    handle to spect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x
global fs
global t
Pw=get(handles.Pw,'String');Pw=str2double(Pw);
Pt=get(handles.Pt,'String');Pt=str2double(Pt);
t0=get(handles.t0,'String');t0=str2double(t0);
T = t0+Pt-Pw-1/fs;  % 接收信号时长
t = 0: 1/fs: T;
L=length(x(1,:));
[s,t,f]=spectrogram(x(1,:),L/10,L*3/40,L/10,fs,'yaxis');
a=[0 T];
b=[0 fs/2];
image(handles.f_wave,a,b,abs(s));%画语谱图
xlabel(handles.f_wave,'时间/s');
ylabel(handles.f_wave,'频率/Hz');
xlim(handles.f_wave,[0,T]);
title(handles.f_wave,'接收信号的时频图');


% --- Executes on button press in estimate.
function estimate_Callback(hObject, eventdata, handles)
% hObject    handle to estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g=estimate;
set(g,'Visible','on');
fig=openfig('estimate.fig');
handles = guihandles(fig);
guidata(fig, handles);


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
newAxes = copyobj(handles.waveform,newFig);   %将axes1中的图复制到新建的figure中
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
newAxes = copyobj(handles.f_wave,newFig);   %将axes1中的图复制到新建的figure中
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
newAxes = copyobj(handles.multi_wave,newFig);   %将axes1中的图复制到新建的figure中
set(newAxes,'Units','default','Position','default');    % 设置图显示的位置
[filename,pathname] = uiputfile({ '*.jpg','figure type(*.jpg)'}, '保存原始波形');
if isequal(filename,0)||isequal(pathname,0)%如果用户选择“取消”，则退出
    return;
else
    fpath=fullfile(pathname,filename);
end
print(newFig,'-djpeg','-r600',fpath);
