function varargout = GUI(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;
lighfield= load('pen.mat');
lighfield =im2double(lighfield.LF(:,:,:,:,1:3));
axes(handles.axes1);
imshow(squeeze(lighfield(6,6,1:375,1:375,1:3)));
set(handles.edit1,'BackgroundColor','green'); 
% Update handles structure
guidata(hObject, handles);
waitforbuttonpress;
set(handles.edit1,'BackgroundColor',[1 1 1]); 
set(handles.pushbutton1,'BackgroundColor','green'); 
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton1,'BackgroundColor',[0.19 0.19 0.19]); 
h = waitbar(0,'please wait...');
lighfield= load('pen.mat');
lighfield =im2double(lighfield.LF(:,:,:,:,1:3));

H =[0.0002         0         0         0   -0.0013
         0    0.0002         0         0   -0.0013
         0         0    0.0004         0   -0.0720
         0         0         0    0.0004   -0.0720
         0         0         0         0    1.0000];
f =     0.1999;

% Extracting SIFT Features For the Centeral Sub-aperture

edge_thresh = 2; % you will have more SIFT features if you use a higher threshold.
IM_center = rgb2gray(squeeze(lighfield(6,6,:,:,:)));

[SIFT_Loc,SIFT_des] = vl_sift(single(IM_center), 'edgethresh', edge_thresh) ; % The outputs will be SIFT Locations and SIFT decriptors 
SIFT_Loc = SIFT_Loc(1:2,:)'; % The First and second rows store the X, Y locations of each SIFT feature respectively. see vl_sift for more details.

% SIFT Matching between SIFT From the Centeral Sub-aperture and all other Sub-apertures


N_SIT = size(SIFT_Loc,1);
match = zeros(N_SIT,25*4);
cnt = 1;
% Note that we use only 25 sub-apertures to find the SIFT matching, outer
% sub-aperture are almost useless due to Vignetting problem.
% and also using lesser sub-aprture will yeid more computational speed.
for i = 4:8
    for j = 4:8
        waitbar((i+j)/16)
        gray = rgb2gray(squeeze(lighfield(j,i,:,:,:)));
        % calculating SIFT for each Sub-apertures
        [SIFTloc_j,SIFTdes_j] = vl_sift(single(gray), 'edgethresh', edge_thresh) ;
        
        % finding the indices of corrspoding SIFT features
        [matchi, matchj] = matchSIFTdesImagesBidirectional(SIFT_des, SIFTdes_j);
        
        % Storing Location of SIFT corrspondences with their sub-aperture indices  
        % each row is corresponding to a specific SIFT feature
        match(matchi,cnt:cnt+3)= [repmat([i j],length(matchi),1) SIFTloc_j(1:2,matchj)' ];
        cnt = cnt+4;
    end
end
% Calculating the 3D Location of each SIFT Features
P3D = zeros(N_SIT,3);
for counter = 1: N_SIT
    
    corres = reshape(match(counter,:),4,25)';
    Nz = find(corres(:,1)~=0); % it might be in some sub-apetures we cannot find correspondencing features. 
    
    % Note that we calculate the 3D location only if the number of correspondences feature are more than 10 
    % less than this number will result in inaccurate 3D Location.
    if length(Nz)>10      
        rays  = H*[corres(Nz,:) ones(length(Nz),1)]'; % transforming pixel coordinates to ray coordinates using intrinsic matrix H
        s = Point_3Dlocation_estimate( rays );
        fz1 =  -s(1)+1;
        P3D(counter,:) = [s(2)/fz1 s(3)/fz1 f/fz1];
    end
end
%
close(h)

K = [SIFT_Loc P3D(:,3)];
X = SIFT_Loc(:,1);
Z = P3D(:,3);
LF1 = lighfield;
relx = [];
relz = [];
d=0;

a = str2num(get(handles.edit1,'String'));

for i = 1:a
    set(handles.text2,'String',['Select object',a])
    rect = getrect(handles.axes1);
    rect = [rect(1) rect(1)+rect(3)];
    temp = find(X<rect(2));
    temp1 = Z(temp);
    temp2 = X(temp);
    temp1 = temp1(find(temp2>rect(1)));
    temp1 = temp1(find(temp1>0));
    temp1 = temp1(find(temp1<0.19*2));
    z = median(temp1);
    relx = [relx rect(1):rect(2)];
    relz = [relz z*ones(1,ceil(rect(2)-rect(1))+1)];
end
h = waitbar(0,'calculating');
% for i = 1:n
%     figure(1)
%     imshow(squeeze(LF1(6,6,:,:,1:3)));
%     title(['select object number',i])
%     rect = getrect;
%     rect = [rect(1) rect(1)+rect(3)];
%     temp = find(X<rect(2));
%     temp1 = Z(temp);
%     temp2 = X(temp);
%     temp1 = temp1(find(temp2>rect(1)));
%     temp1 = temp1(find(temp1>0));
%     temp1 = temp1(find(temp1<0.19*2));
%     z = median(temp1);
%     relx = [relx rect(1):rect(2)];
%     relz = [relz z*ones(1,ceil(rect(2)-rect(1))+1)];
% end
[~,~,x,y,~] = size(LF1);
[relx,I]=sort(relx,'ascend');
relz = relz(I);
mini = min(relx);
maxi = max(relx);
relx = [1:mini-1 relx maxi+1:x];
relz = [relz(1)*ones(1,ceil(mini-1)) relz relz(end)*ones(1,ceil(x-maxi))];
H1 = [0.0004 0 -0.072 ; 0 0.0004 -0.072 ; 0 0 1];
A = H1*[relx ;ones(size(relx));ones(size(relx))];
relu = A(1,:);
%%
H11 = inv(H1);
H2 = [0.0002 0 -0.0013 ; 0 0.0002 -0.0013 ; 0 0 1];
%%%%%
[l,k] = ndgrid(1:x,1:y);
k = reshape(k,[1,x*y]);
l = reshape(l,[1,x*y]);
%%%%%
[i,j] = ndgrid(1:11,1:11);
i = reshape(i,[1,11^2]);
j = reshape(j,[1,11^2]);
%%%%%
A = H1*[k;l;ones(1,x*y)];
u = A(1,:);
v = A(2,:);
%%%%%
A = H2*[i;j;ones(1,121)];
s = A(1,:);
t = A(2,:);
%%%%%
up = u;
vp = v;
%%
temp = reshape(u,[x y]);
temp = temp(1,:);
ztemp = interpn(relu,relz,temp,'spline',0);
%%
zp = reshape(ones(x,1)*ztemp,[1,x*y]);
F=0.19;
LF1 = im2double(LF1);
refocused=zeros(x,y,3);
    for i = 3:9
        for j = 3:9
            for k=1:3
                waitbar((i+j+k)/21)
                CurSlice=squeeze(LF1(j,i,:,:,k));
                u2 = up+s(i)*(1-F./zp);
                v2 = vp+s(j)*(1-F./zp);
                A = H11*[u2;v2;ones(1,x*y)];
                k2 = reshape(A(1,:),[x,y]);
                l2 = reshape(A(2,:),[x,y]);
                CurSlice = interpn(CurSlice,l2,k2,'spline',0);
                refocused(:,:,k)=CurSlice+refocused(:,:,k);
            end
        end
    end   
refocused = refocused/49;
close(h)
axes(handles.axes1);
set(handles.text2,'String','final result')


guidata(hObject, handles);