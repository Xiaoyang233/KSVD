function varargout = KSVD_Show(varargin)
% KSVD_SHOW MATLAB code for KSVD_Show.fig
%      KSVD_SHOW, by itself, creates a new KSVD_SHOW or raises the existing
%      singleton*.
%
%      H = KSVD_SHOW returns the handle to a new KSVD_SHOW or the handle to
%      the existing singleton*.
%
%      KSVD_SHOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KSVD_SHOW.M with the given input arguments.
%
%      KSVD_SHOW('Property','Value',...) creates a new KSVD_SHOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KSVD_Show_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KSVD_Show_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KSVD_Show

% Last Modified by GUIDE v2.5 29-Apr-2015 19:56:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KSVD_Show_OpeningFcn, ...
                   'gui_OutputFcn',  @KSVD_Show_OutputFcn, ...
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


% --- Executes just before KSVD_Show is made visible.
function KSVD_Show_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to KSVD_Show (see VARARGIN)

% Choose default command line output for KSVD_Show
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes KSVD_Show wait for user response (see UIRESUME)
% uiwait(handles.figure1);

ImageFile  = 'opening1.jpg';
ScreenSize = get(0,'ScreenSize');
jImage     = im2java(imread(ImageFile));
jfBounds   = num2cell([...
    (ScreenSize(3)-jImage.getWidth)/2 ...
    (ScreenSize(4)-jImage.getHeight)/2 ...
    jImage.getWidth ...
    jImage.getHeight]);
jFrame     = javax.swing.JFrame;
icon       = javax.swing.ImageIcon(jImage);
label      = javax.swing.JLabel(icon);
jFrame.getContentPane.add(label);
jFrame.setUndecorated(true)
jFrame.setBounds(jfBounds{:});
jFrame.pack
jFrame.show
pause(5)
jFrame.dispose



% --- Outputs from this function are returned to the command line.
function varargout = KSVD_Show_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pressRate1 = get(handles.radiobutton2,'value');
pressRate2 = get(handles.radiobutton3,'value');
if pressRate1 == 1
    load Dic_seterr4;
else
    load Dic;
end

testDic();

axes(handles.axes1);
imshow(testPic,[]);
axes(handles.axes2);
imshow(reconstKSVD,[]);
axes(handles.axes4);
 I1 = displayDictionaryElementsAsImage(trainDic, floor(sqrt(K)), floor(size(trainDic,2)/floor(sqrt(K))),8,8,0);
 
 axes(handles.axes3);
 imshow(reconstDCT,[]);
  axes(handles.axes5);
  DctDic=output.D;
 I2 = displayDictionaryElementsAsImage(DctDic, floor(sqrt(K)), floor(size(DctDic,2)/floor(sqrt(K))),8,8,0);
 
set(handles.edit1,'String',['PSNR: ',num2str(PSNROutKsvd),' db']);
set(handles.edit2,'String',['PSNR: ',num2str(PSNROutDct),' db']);

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
%guidata(hObject, handles);

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

ha=axes('units','normalized','position',[0 0 1 1]);
uistack(ha,'down')
II=imread('ppt1.jpg');
image(II)
colormap gray
set(ha,'handlevisibility','off','visible','off');





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


% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 I=imread('button.jpg');
 h=findobj('tag','pushbutton1');
set(h,'CData', I);
