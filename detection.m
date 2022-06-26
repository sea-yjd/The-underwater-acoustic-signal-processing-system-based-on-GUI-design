function varargout = detection(varargin)
% DETECTION MATLAB code for detection.fig
%      DETECTION, by itself, creates a new DETECTION or raises the existing
%      singleton*.
%
%      H = DETECTION returns the handle to a new DETECTION or the handle to
%      the existing singleton*.
%
%      DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTION.M with the given input arguments.
%
%      DETECTION('Property','Value',...) creates a new DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detection

% Last Modified by GUIDE v2.5 19-Jun-2022 15:40:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @detection_OpeningFcn, ...
    'gui_OutputFcn',  @detection_OutputFcn, ...
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


% --- Executes just before detection is made visible.
function detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to detection (see VARARGIN)

% Choose default command line output for detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = detection_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Once.
function Once_Callback(hObject, eventdata, handles)
% hObject    handle to Once (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Once


% --- Executes on button press in Many.
function Many_Callback(hObject, eventdata, handles)
% hObject    handle to Many (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Many



function time_Callback(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time as text
%        str2double(get(hObject,'String')) returns contents of time as a double


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fs,x,t,s,f0,Bw,Pw,t0,theata,N,s_s,LFM,CW,SNR_db]=part1;
Once=get(handles.Once,'Value');
Many=get(handles.Many,'Value');
M=Bw/Pw;
Pf=0.1;
Nk=1000;
%实验曲线
if Once==1
    n_x = [1:Nk];
    A = sqrt(0.5); %令信号的功率始终为一
    noisepower = 10^(-SNR_db/10);
    Es = Nk;
    T = sqrt(noisepower*Es)*norminv(1-Pf/2);
        if LFM==1
            sig=A*cos(2*pi*((f0-Bw/2-M)*n_x+0.5*M*n_x.^2));
        else
            sig=A*1/sqrt(Pw)*cos(2*pi*f0*n_x);
        end
%         st = A*cos(2*pi*f0.*n_x);
        xt = awgn(sig,SNR_db,'measured');
        %             st = 1/sqrt(Nk)*rectpuls(n_x,Nk).*exp(j*2*pi*f0*n_x);
        %         xt=st+noise;
        %         Qx=sqrt((sin(2*pi*f0*n_x)*(xt)')^2+(cos(2*pi*f0*n_x)*(xt)')^2);
        G = sig*xt';
        if G>T
           set(handles.det,'String','是');
        else
           set(handles.det,'String','否');
        end
else
NS=get(handles.time,'String');NS=str2double(NS);
SNR=[-35:0.5:15];

% for e=1:length(Nk)
Pd=[];

for f=1:length(SNR)
    Nm=0;
    %         phi=rand(1)*2*pi;
    r = SNR(f);
    %         A=sqrt(2*10^(r/10)); %方差为1
    A = sqrt(0.5); %令信号的功率始终为一
    noisepower = 10^(-r/10);
    Es = Nk;
    T = sqrt(noisepower*Es)*norminv(1-Pf/2);
    for g=1:NS
        %         phi=rand(1)*2*pi;
        %         noise=randn(1,Nk);
        n_x = [1:Nk];
        if LFM==1
            sig=A*cos(2*pi*((f0-Bw/2-M)*n_x+0.5*M*n_x.^2));
        else
            sig=A*1/sqrt(Pw)*cos(2*pi*f0*n_x);
        end
%         st = A*cos(2*pi*f0.*n_x);
        xt = awgn(sig,r,'measured');
        %             st = 1/sqrt(Nk)*rectpuls(n_x,Nk).*exp(j*2*pi*f0*n_x);
        %         xt=st+noise;
        %         Qx=sqrt((sin(2*pi*f0*n_x)*(xt)')^2+(cos(2*pi*f0*n_x)*(xt)')^2);
        G = sig*xt';
        if G>T
            Nm=Nm+1;
        end
    end
    Pd(f)=Nm/NS;
end
set(handles.pd,'String',num2str(Pd(1,91)));
plot(handles.axes1,SNR,Pd,'LineWidth',2);hold on
xlabel(handles.axes1,'输入信噪比SNR/dB');ylabel(handles.axes1,'检测概率Pd');
title(handles.axes1,'检测概率曲线(改变样本数)');
end
% if LFM==1
%     sig=cos(2*pi*((f0-Bw/2-M)*t+0.5*M*t.^2));
% else
%     sig=1/sqrt(Pw)*cos(2*pi*f0*t);
% end
% if Once==1
%     noi = x(1,:)-sig;
%     noisevar = noi*noi'/2;
%     [T,G] = signaldetection(x(1,:),sig,noisevar);
%     if G>double(T)
%         set(handles.det,'String','是');
%     else
%         set(handles.det,'String','否');
%     end
% else
%     Nmot=get(handles.time,'String');Nmot=str2double(Nmot);
%     SNR=-10:0.5:15;
%     pf = 0.0001;
%     Es = sig*sig'; 
%     
%     for k=1:length(SNR)
%         %function [therhold,G]=signaldetection(x,s,omega)
%         x(1,:)=awgn(s(1,:),SNR(k),'measured'); 
%         noi = x(1,:)-s(1,:);
%         noisevar = noi*noi'/2;
%         t=0;
%         syms eta;
%         eta = solve(2*(1-normcdf(eta*sqrt(1/noisevar/Es)))==pf,eta);
%         T = eta;
%         for i=1:Nmot
%             x(i,:)=awgn(s(1,:),SNR(k),'measured'); 
%             G = x(i,:)*sig';
%            
%             disp(['SNR:',num2str(SNR(k)),'次数:',num2str(i),':::',num2str(t)]);
%             if G>double(T)
%                 t=t+1;
% %             else
% %                 t=t;
%             end
%         end
%         Pd(1,k)=t/Nmot;
%     end
%     plot(handles.axes1,SNR,Pd,'LineWidth',2);hold on
%     xlabel(handles.axes1,'输入信噪比SNR（dB））');ylabel(handles.axes1,'检测概率Pd');
%     title('检测概率曲线');
%     %legend('理论','实验');
%     set(handles.pd,'String',num2str(Pd(1,3)));
% end

% H1=1/(2*sqrt(pi))*exp(-(x-s_s).^2/2);
% H0=1/(2*sqrt(pi))*exp(-x.^2/2);
% lambda_x=H1/H0;
% Pf=0.1;
% if Many==1
% Ntrial=get(handles.time,'String');Ntrial=str2double(Ntrial);
% mf = 1;
% y = mf'*x(1,:);  % apply the matched filter
% z = real(y);
% Pfa = 1e-3;
% snrthreshold = db2pow(npwgnthresh(Pfa, 1,'coherent'));
% mfgain = mf'*mf;
% % To match the equation in the text above
% % npower - N
% % mfgain - M
% % snrthreshold - SNR
% spower=s_s*s_s';
% SNR = db2pow(SNR_db);      % SNR in linear scale
% npower=spower/SNR;
% threshold = sqrt(npower*mfgain*snrthreshold);
% Pd = sum(z*z'>threshold)/Ntrial
% namp = sqrt(npower/2);
% n = namp.*randn(1,length(x))+1i.*randn(1,length(x));
% y = mf'*n;
% z = real(y);
% Pfa = sum(z*z'>threshold)/Ntrial
% %axes1
% rocsnr(SNR_db,'SignalType','NonfluctuatingCoherent','MinPfa',1e-4);
% set(handles.Pd,'String',num2str(Pd));
% set(handles.Pfa,'String',num2str(Pfa));
% else
%     mf = 1;
% y = mf'*x(1,:);  % apply the matched filter
% z = real(y);
% Pfa = 1e-3;
% snrthreshold = db2pow(npwgnthresh(Pfa, 1,'coherent'));
% mfgain = mf'*mf;
% % To match the equation in the text above
% % npower - N
% % mfgain - M
% % snrthreshold - SNR
% spower=s_s*s_s';
% SNR = db2pow(SNR_db);      % SNR in linear scale
% npower=spower/SNR;
% threshold = sqrt(npower*mfgain*snrthreshold);
% Pd = sum(z*z'>threshold)/Ntrial;
% if Pd>=0
%     set(handles.det,'String','是');
% else
%     set(handles.det,'String','否');
% end
% end
% --- Executes on button press in detection_b.
function detection_b_Callback(hObject, eventdata, handles)
% hObject    handle to detection_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fs,x,t,s,f0,Bw,Pw,t0,theata,N,s_s,LFM,CW,SNR_db,B]=part1;
Once_b=get(handles.Once_b,'Value');
Many_b=get(handles.Many_b,'Value');
M=Bw/Pw;
Pf=0.1;
Nk=1000;
%实验曲线
if Once_b==1
    n_x = [1:Nk];
    A = sqrt(0.5); %令信号的功率始终为一
    noisepower = 10^(-SNR_db/10);
    Es = Nk;
    T = sqrt(noisepower*Es)*norminv(1-Pf/2);
    %信号
    % 阵元上的其他噪声
    noise_env = randn(N,length(n_x));
    
    s1_amplitude = A*randn(1);  % 线谱幅度
    s1_phase = rand(1) * 2 * pi;  % 线谱相位
    s1_t = s1_amplitude * cos(2*pi*f0*n_x+s1_phase);  % 线谱成分
    % 宽带噪声
    fl = f0 - B/2;
    fh = f0 + B/2;
    s2_t = randn(1,length(n_x));
    s2_fft = fft(s2_t);
    num_fl = floor(fl/fs*length(n_x));
    num_fh = ceil(fh/fs*length(n_x));
    s2_fft = [zeros(1,num_fl-1) s2_fft(num_fl:num_fh) zeros(1,length(n_x)-num_fh)];
    s2_t = real(ifft(s2_fft));
    %海洋环境噪声
    s_t = s1_t + s2_t;
    power_s = 10^(SNR_db/10) * var(noise_env(1,:));
    sig = sqrt(power_s) / std(s_t) * s_t;
%         st = A*cos(2*pi*f0.*n_x);
        xt = awgn(sig,SNR_db,'measured');
        %             st = 1/sqrt(Nk)*rectpuls(n_x,Nk).*exp(j*2*pi*f0*n_x);
        %         xt=st+noise;
        %         Qx=sqrt((sin(2*pi*f0*n_x)*(xt)')^2+(cos(2*pi*f0*n_x)*(xt)')^2);
        G = sig*xt';
        if G>T
           set(handles.det_b,'String','是');
        else
           set(handles.det_b,'String','否');
        end
else
NS=get(handles.time_b,'String');NS=str2double(NS);
SNR=[-35:0.5:15];

% for e=1:length(Nk)
Pd=[];

for f=1:length(SNR)
    Nm=0;
    %         phi=rand(1)*2*pi;
    r = SNR(f);
    %         A=sqrt(2*10^(r/10)); %方差为1
    A = sqrt(0.5); %令信号的功率始终为一
    noisepower = 10^(-r/10);
    Es = Nk;
    T = sqrt(noisepower*Es)*norminv(1-Pf/2);
    for g=1:NS
        %         phi=rand(1)*2*pi;
        %         noise=randn(1,Nk);
        n_x = [1:Nk];
            %信号
            % 阵元上的其他噪声
            noise_env = randn(N,length(n_x));
            
            s1_amplitude = A*randn(1);  % 线谱幅度
            s1_phase = rand(1) * 2 * pi;  % 线谱相位
            s1_t = s1_amplitude * cos(2*pi*f0*n_x+s1_phase);  % 线谱成分
            % 宽带噪声
            fl = f0 - B/2;
            fh = f0 + B/2;
            s2_t = randn(1,length(n_x));
            s2_fft = fft(s2_t);
            num_fl = floor(fl/fs*length(n_x));
            num_fh = ceil(fh/fs*length(n_x));
            s2_fft = [zeros(1,num_fl-1) s2_fft(num_fl:num_fh) zeros(1,length(n_x)-num_fh)];
            s2_t = real(ifft(s2_fft));
            %海洋环境噪声
            s_t = s1_t + s2_t;
            power_s = 10^(SNR_db/10) * var(noise_env(1,:));
            sig = sqrt(power_s) / std(s_t) * s_t;
%         st = A*cos(2*pi*f0.*n_x);
        xt = awgn(sig,r,'measured');
        %             st = 1/sqrt(Nk)*rectpuls(n_x,Nk).*exp(j*2*pi*f0*n_x);
        %         xt=st+noise;
        %         Qx=sqrt((sin(2*pi*f0*n_x)*(xt)')^2+(cos(2*pi*f0*n_x)*(xt)')^2);
        G = sig*xt';
        if G>T
            Nm=Nm+1;
        end
    end
    Pd(f)=Nm/NS;
end
set(handles.pd_b,'String',num2str(Pd(1,91)));
plot(handles.axes2,SNR,Pd,'LineWidth',2);hold on
xlabel(handles.axes2,'输入信噪比SNR/dB');ylabel(handles.axes1,'检测概率Pd');
title(handles.axes2,'检测概率曲线(改变样本数)');
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); %指定需要清空的坐标轴
cla reset;
axes(handles.axes2); %指定需要清空的坐标轴
cla reset;
handles.waveform.Toolbar.Visible = 'on';




function det_Callback(hObject, eventdata, handles)
% hObject    handle to det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of det as text
%        str2double(get(hObject,'String')) returns contents of det as a double


% --- Executes during object creation, after setting all properties.
function det_CreateFcn(hObject, eventdata, handles)
% hObject    handle to det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd_Callback(hObject, eventdata, handles)
% hObject    handle to pd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd as text
%        str2double(get(hObject,'String')) returns contents of pd as a double


% --- Executes during object creation, after setting all properties.
function pd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_b_Callback(hObject, eventdata, handles)
% hObject    handle to time_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_b as text
%        str2double(get(hObject,'String')) returns contents of time_b as a double


% --- Executes during object creation, after setting all properties.
function time_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function det_b_Callback(hObject, eventdata, handles)
% hObject    handle to det_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of det_b as text
%        str2double(get(hObject,'String')) returns contents of det_b as a double


% --- Executes during object creation, after setting all properties.
function det_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to det_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pd_b_Callback(hObject, eventdata, handles)
% hObject    handle to pd_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pd_b as text
%        str2double(get(hObject,'String')) returns contents of pd_b as a double


% --- Executes during object creation, after setting all properties.
function pd_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pd_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
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
