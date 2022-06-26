function varargout = estimate(varargin)
% ESTIMATE MATLAB code for estimate.fig
%      ESTIMATE, by itself, creates a new ESTIMATE or raises the existing
%      singleton*.
%
%      H = ESTIMATE returns the handle to a new ESTIMATE or the handle to
%      the existing singleton*.
%
%      ESTIMATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ESTIMATE.M with the given input arguments.
%
%      ESTIMATE('Property','Value',...) creates a new ESTIMATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before estimate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to estimate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help estimate

% Last Modified by GUIDE v2.5 10-Jun-2022 20:56:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @estimate_OpeningFcn, ...
                   'gui_OutputFcn',  @estimate_OutputFcn, ...
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


% --- Executes just before estimate is made visible.
function estimate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to estimate (see VARARGIN)

% Choose default command line output for estimate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes estimate wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = estimate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1); %指定需要清空的坐标轴
cla reset;
axes(handles.axes2); %指定需要清空的坐标轴
cla reset;


function Dis_v_Callback(hObject, eventdata, handles)
% hObject    handle to Dis_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dis_v as text
%        str2double(get(hObject,'String')) returns contents of Dis_v as a double


% --- Executes during object creation, after setting all properties.
function Dis_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dis_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DoA_v_Callback(hObject, eventdata, handles)
% hObject    handle to DoA_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DoA_v as text
%        str2double(get(hObject,'String')) returns contents of DoA_v as a double


% --- Executes during object creation, after setting all properties.
function DoA_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DoA_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_o_v_Callback(hObject, eventdata, handles)
% hObject    handle to v_o_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v_o_v as text
%        str2double(get(hObject,'String')) returns contents of v_o_v as a double


% --- Executes during object creation, after setting all properties.
function v_o_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_o_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fs,x,t,s,f0,Bw,Pw,t0,theata,N,s_s,LFM,CW,SNR]=part1;
Once=get(handles.Once_b,'Value');
Many=get(handles.Many_b,'Value');
c = 1500;% 声速（m/s）
lambda0 = c / f0;% 波长
d = lambda0 / 2;% 阵元间距为信号波长的一半（标准线列阵 SLA）
theataTArray = [0:1:180-1]';
uTArray = cos(theataTArray * pi / 180);
theataTLen = length(theataTArray);
% 扫描角（角度空间和u空间）
n=1:N;
w=ones(N,1)/N;
% %% 阵元上的信号
% delaynum_object = round((0:N-1).'*d*cos(theata)/c*fs);  % 以第一个阵元上的信号为参考信号   %% 
% v_object = exp(-1j*2*pi/length(t)*delaynum_object*(0:length(t)-1));  
% s_fft = fft(s);
% s_fft_vector = v_object .* s_fft;
% s_t_vector = ifft(s_fft_vector.').';
% x_t_vector = s_t_vector + randn(N,length(t));
% x_fft_vector = fft(x_t_vector.').';
%单次估计
if Once==1
    
for theataTSelect = 1 : theataTLen
% 波束图

theataT = theataTArray(theataTSelect, 1);
phy_scan=2*d*pi*((cosd(theata)-cosd(theataTArray(theataTSelect))))/lambda0;
A_scan=exp(1i*(((n-1).'-(N-1)/2))*phy_scan');
%A_scan=exp(1i*((n-1).')*phy_scan);      %阵列流形
noise=sqrt(1)*randn(1,length(s(1,:)));    
wn=w.*A_scan;
%x=sum(wn*fft(s_s));%
%s_fft=fft(s_s);
%x_x=x(1,:)+noise;
y(theataTSelect,:)=sum(conj(wn).*fft(x(1,:)),1);%+noise
yifft(theataTSelect,:)=ifft(y(theataTSelect,:));%输出再加上噪声
p(theataTSelect)=yifft(theataTSelect,:)*yifft(theataTSelect,:)';
%y(thetaTSelect,:)=ifft(Y(thetaTSelect,:));
end

%功率

[~,D]=max(p);
DoA=D-1;
%DoA估计
% plot(handles.axes1,t,x_x,'linewidth',2);
% xlabel(handles.axes1,'接收信号的时间/s');
% ylabel(handles.axes1,'信号幅值');
% title(handles.axes1,'阵元1接收的被动信号波形');
% plot(handles.axes2,t,real(yifft(DoA,:)),'linewidth',1.5);
% xlabel(handles.axes2,'接收信号的时间/s');
% ylabel(handles.axes2,'信号幅值');
% title(handles.axes2,'均匀加权后输出波形');
% set(handles.DoA_b,'String',num2str(DoA));
%均匀加权后输出波形

 
 
%%多次统计
else
    times=get(handles.times,'String');times=str2double(times);
    for i=1:times
        for theataTSelect = 1 : theataTLen
% 波束图
        theataT = theataTArray(theataTSelect, 1);
        phy_scan=2*d*pi*((cosd(theata)-cosd(theataTArray(theataTSelect))))/lambda0;
        A_scan=exp(1i*(((n-1).'-(N-1)/2))*phy_scan');
%A_scan=exp(1i*((n-1).')*phy_scan);      %阵列流形
%noise=sqrt(1)*randn(1,length(s(1,:)));    
        wn=w.*A_scan;
%x=sum(wn*fft(s_s));%
%s_fft=fft(s_s);
        x(1,:)=awgn(s(1,:),SNR,'measured');
        y(theataTSelect,:)=sum(conj(wn).*fft(x(1,:)),1);%
        yifft(theataTSelect,:)=ifft(y(theataTSelect,:));%输出再加上噪声
        p(theataTSelect)=yifft(theataTSelect,:)*yifft(theataTSelect,:)';
%y(thetaTSelect,:)=ifft(Y(thetaTSelect,:));
        end

%功率
%DoA估计
[~,D]=max(p);
DoA(1,i)=D-1;

    end
%%均值与方差
DoA_u=mean(DoA);
DoA_v=var(DoA,0);
set(handles.DoA_u_b,'String',num2str(DoA_u));
set(handles.DoA_v_b,'String',num2str(DoA_v));

end

% --- Executes on button press in Once_b.
function Once_b_Callback(hObject, eventdata, handles)
% hObject    handle to Once_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Once_b


% --- Executes on button press in Many_b.
function Many_b_Callback(hObject, eventdata, handles)
% hObject    handle to Many_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Many_b



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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



function times_Callback(hObject, eventdata, handles)
% hObject    handle to times (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of times as text
%        str2double(get(hObject,'String')) returns contents of times as a double


% --- Executes during object creation, after setting all properties.
function times_CreateFcn(hObject, eventdata, handles)
% hObject    handle to times (see GCBO)
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
[fs,x,t,s,f0,Bw,Pw,t0,theata,N,s_s,LFM,CW,SNR]=part1;
Once=get(handles.Once,'Value');
Many=get(handles.Many,'Value');
c = 1500;% 声速（m/s）
lambda0 = c / f0;% 波长
d = lambda0 / 2;% 阵元间距为信号波长的一半（标准线列阵 SLA）
theata0 = 1:180;%弧度换算角度
theataTArray = [0:1:180-1]';
uTArray = cos(theataTArray * pi / 180);
theataTLen = length(theataTArray);
% 扫描角（角度空间和u空间）
n=1:N;
w=ones(N,1)/N;
%单次估计
if Once==1
    
for theataTSelect = 1 : theataTLen
% 波束图

theataT = theataTArray(theataTSelect, 1);
phy_scan=2*d*pi*((cosd(theata)-cosd(theataTArray(theataTSelect))))/lambda0;
A_scan=exp(1i*(((n-1).'-(N-1)/2))*phy_scan');
%A_scan=exp(1i*((n-1).')*phy_scan);      %阵列流形
%noise=sqrt(1)*randn(1,length(s(1,:)));    
wn=w.*A_scan;
%x=sum(wn*fft(s_s));%
%s_fft=fft(s_s);
y(theataTSelect,:)=sum(conj(wn).*fft(x(1,:)),1);%
yifft(theataTSelect,:)=ifft(y(theataTSelect,:));%输出再加上噪声
p(theataTSelect)=yifft(theataTSelect,:)*yifft(theataTSelect,:)';
%y(thetaTSelect,:)=ifft(Y(thetaTSelect,:));
end

%功率


[~,D]=max(p);
DoA=D-1;
% plot(handles.axes1,t,x(1,:),'linewidth',2);
% xlabel(handles.axes1,'接收信号的时间/s');
% ylabel(handles.axes1,'信号幅值');
% title(handles.axes1,'阵元1接收的被动信号波形');
% plot(handles.axes2,t,real(yifft(DoA,:)),'linewidth',1.5);
% xlabel(handles.axes2,'接收信号的时间/s');
% ylabel(handles.axes2,'信号幅值');
% title(handles.axes2,'均匀加权后输出波形');
%DoA估计
set(handles.DoA,'String',num2str(DoA));
%均匀加权后输出波形
 %发射信号
 t1=0:1/fs:Pw-1/fs;
 %t2=Pw:1/fs:1.4-1/fs;
 %t_t=[t1,t2];
 %s_o=[s_s,zeros(1,length(t2))];
 %plot(handles.axes3,t1,s_s,'linewidth',1.5);
 [a1,b1]=xcorr(yifft(DoA,:),s_s);
 [~,I1] = max(abs(a1));
lagDiff1 = b1(I1);
timeDiff1 = lagDiff1/fs;
Dis=timeDiff1*c/2;
%目标距离
set(handles.Dis,'String',num2str(Dis));
%目标速度
f_num=10;
M=Bw/Pw;

t_tao=t1+timeDiff1;
for v=1:f_num
    %sig_LFM=cos(2*pi*((f0-Bw/2-M)*t_tao+0.5*M*t_tao.^2));
    if LFM==1
    df=0.33*Bw;
    sig_LFM_scan=cos(2*pi*((f0+(v-1)*df-Bw/2-M)*t_tao+0.5*M*t_tao.^2));
    sig_scan_n=awgn(sig_LFM_scan,SNR,'measured');
    else
    df=0.44/Pw;
    sig_CW_scan=1/sqrt(Pw)*cos(2*pi*(f0+(v-1)*df)*t_tao);
    sig_scan_n=awgn(sig_CW_scan,SNR,'measured');
    end
    [a2,b2]=xcorr(sig_scan_n,s_s);
    [I2,~] = max(abs(a2));
    f_max(1,v)=I2;
end
[~,v_index]=max(f_max);
v_o=lambda0*(v_index-1)*df;
set(handles.v_o,'String',num2str(v_o));


%%多次统计
else
    times=get(handles.times,'String');times=str2double(times);
    for i=1:times
        for theataTSelect = 1 : theataTLen
% 波束图
        theataT = theataTArray(theataTSelect, 1);
        phy_scan=2*d*pi*((cosd(theata)-cosd(theataTArray(theataTSelect))))/lambda0;
        A_scan=exp(1i*(((n-1).'-(N-1)/2))*phy_scan');
%A_scan=exp(1i*((n-1).')*phy_scan);      %阵列流形
%noise=sqrt(1)*randn(1,length(s(1,:)));    
        wn=w.*A_scan;
%x=sum(wn*fft(s_s));%
%s_fft=fft(s_s);
        x(1,:)=awgn(s(1,:),SNR,'measured');
        y(theataTSelect,:)=sum(conj(wn).*fft(x(1,:)),1);%
        yifft(theataTSelect,:)=ifft(y(theataTSelect,:));%输出再加上噪声
        p(theataTSelect)=yifft(theataTSelect,:)*yifft(theataTSelect,:)';
%y(thetaTSelect,:)=ifft(Y(thetaTSelect,:));
        end

%功率
%DoA估计
[~,D]=max(p);
DoA(1,i)=D-1;
%均匀加权后输出波形
 %发射信号
 t1=0:1/fs:Pw-1/fs;
 %t2=Pw:1/fs:1.4-1/fs;
 %t_t=[t1,t2];
 %s_o=[s_s,zeros(1,length(t2))];
 [a1,b1]=xcorr(yifft(DoA(1,i),:),s_s);
 [~,I1] = max(abs(a1));
lagDiff1 = b1(I1);
timeDiff1 = lagDiff1/fs;
%目标距离;
Dis(1,i)=timeDiff1*c/2;

%目标速度
f_num=10;
M=Bw/Pw;

t_tao=t1+timeDiff1;
for v=1:f_num
    %sig_LFM=cos(2*pi*((f0-Bw/2-M)*t_tao+0.5*M*t_tao.^2));
    if LFM==1
    df=0.33*Bw;
    sig_LFM_scan=cos(2*pi*((f0+(v-1)*df-Bw/2-M)*t_tao+0.5*M*t_tao.^2));
    sig_scan_n=awgn(sig_LFM_scan,SNR,'measured');
    else
    df=0.44/Pw;
    sig_CW_scan=1/sqrt(Pw)*cos(2*pi*(f0+(v-1)*df)*t_tao);
    sig_scan_n=awgn(sig_CW_scan,SNR,'measured');
    end
    [a2,b2]=xcorr(sig_scan_n,s_s);
    [I2,~] = max(abs(a2));
    f_max(1,v)=I2;
end
[~,v_index]=max(f_max);
v_o(1,i)=lambda0*(v_index-1)*df;
    end
%%均值与方差
Dis_u=mean(Dis);
v_o_u=mean(v_o);
DoA_u=mean(DoA);
Dis_v=var(Dis,0);
v_o_v=var(v_o,0);
DoA_v=var(DoA,0);
set(handles.Dis_u,'String',num2str(Dis_u));
set(handles.v_o_u,'String',num2str(v_o_u));
set(handles.DoA_u,'String',num2str(DoA_u));
set(handles.Dis_v,'String',num2str(Dis_v));
set(handles.v_o_v,'String',num2str(v_o_v));
set(handles.DoA_v,'String',num2str(DoA_v));

end
    

 
%     xlabel('t(s)','fontsize',10);
%     ylabel('幅值')
%     xlim([0,2]); ylim([-3,3]);
%     title(['扫描角\theta_T = ', num2str(thetaTArray(k)), '\circ'],'fontsize',1);
%     set(gca,'xtick',0:0.5:2,'ytick',-3:1:3)
%     set(gca,'linewidth',1,'fontsize',10);
%多次估计
% else if Many==1
 %end


function Dis_Callback(hObject, eventdata, handles)
% hObject    handle to Dis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dis as text
%        str2double(get(hObject,'String')) returns contents of Dis as a double


% --- Executes during object creation, after setting all properties.
function Dis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DoA_Callback(hObject, eventdata, handles)
% hObject    handle to DoA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DoA as text
%        str2double(get(hObject,'String')) returns contents of DoA as a double


% --- Executes during object creation, after setting all properties.
function DoA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DoA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_o_Callback(hObject, eventdata, handles)
% hObject    handle to v_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v_o as text
%        str2double(get(hObject,'String')) returns contents of v_o as a double


% --- Executes during object creation, after setting all properties.
function v_o_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Dis_u_Callback(hObject, eventdata, handles)
% hObject    handle to Dis_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dis_u as text
%        str2double(get(hObject,'String')) returns contents of Dis_u as a double


% --- Executes during object creation, after setting all properties.
function Dis_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dis_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DoA_u_Callback(hObject, eventdata, handles)
% hObject    handle to DoA_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DoA_u as text
%        str2double(get(hObject,'String')) returns contents of DoA_u as a double


% --- Executes during object creation, after setting all properties.
function DoA_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DoA_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_o_u_Callback(hObject, eventdata, handles)
% hObject    handle to v_o_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v_o_u as text
%        str2double(get(hObject,'String')) returns contents of v_o_u as a double


% --- Executes during object creation, after setting all properties.
function v_o_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_o_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DoA_b_Callback(hObject, eventdata, handles)
% hObject    handle to DoA_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DoA_b as text
%        str2double(get(hObject,'String')) returns contents of DoA_b as a double


% --- Executes during object creation, after setting all properties.
function DoA_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DoA_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DoA_u_b_Callback(hObject, eventdata, handles)
% hObject    handle to DoA_u_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DoA_u_b as text
%        str2double(get(hObject,'String')) returns contents of DoA_u_b as a double


% --- Executes during object creation, after setting all properties.
function DoA_u_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DoA_u_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DoA_v_b_Callback(hObject, eventdata, handles)
% hObject    handle to DoA_v_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DoA_v_b as text
%        str2double(get(hObject,'String')) returns contents of DoA_v_b as a double


% --- Executes during object creation, after setting all properties.
function DoA_v_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DoA_v_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
