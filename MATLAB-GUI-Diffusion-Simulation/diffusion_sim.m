function varargout = diffusion_sim(varargin)
%DIFFUSION_SIM M-file for diffusion_sim.fig
%      DIFFUSION_SIM, by itself, creates a new DIFFUSION_SIM or raises the existing
%      singleton*.
%
%      H = DIFFUSION_SIM returns the handle to a new DIFFUSION_SIM or the handle to
%      the existing singleton*.
%
%      DIFFUSION_SIM('Property','Value',...) creates a new DIFFUSION_SIM using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to diffusion_sim_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DIFFUSION_SIM('CALLBACK') and DIFFUSION_SIM('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DIFFUSION_SIM.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help diffusion_sim

% Last Modified by GUIDE v2.5 08-Sep-2014 16:55:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @diffusion_sim_OpeningFcn, ...
                   'gui_OutputFcn',  @diffusion_sim_OutputFcn, ...
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


% --- Executes just before diffusion_sim is made visible.
function diffusion_sim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for diffusion_sim
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes diffusion_sim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = diffusion_sim_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function initC_Callback(hObject, eventdata, handles)
% hObject    handle to initC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initC as text
%        str2double(get(hObject,'String')) returns contents of initC as a double


% --- Executes during object creation, after setting all properties.
function initC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diff_coef_Callback(hObject, eventdata, handles)
% hObject    handle to diff_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diff_coef as text
%        str2double(get(hObject,'String')) returns contents of diff_coef as a double


% --- Executes during object creation, after setting all properties.
function diff_coef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diff_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pause_tog.
function pause_tog_Callback(hObject, eventdata, handles)
% hObject    handle to pause_tog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pause_tog

  if get(hObject,'Value') == get(hObject,'Max')
      set(hObject,'String','Resume');
      set(handles.info_text,'String','Choose the figure to export as an image');
      set(handles.exp_img,'Enable','on');
      set(handles.exp_fig1,'Enable','on');
      set(handles.exp_fig2,'Enable','on');
      set(handles.exp_fig3,'Enable','on');
      set(handles.exp_gui,'Enable','on');
  else
      set(hObject,'String','Pause');
      set(handles.info_text,'String','');
      set(handles.exp_img,'Enable','off');
      set(handles.exp_fig1,'Enable','off');
      set(handles.exp_fig2,'Enable','off');
      set(handles.exp_fig3,'Enable','off');
      set(handles.exp_gui,'Enable','off');
  end

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;

% --- Executes on button press in exp_vid.
function exp_vid_Callback(hObject, eventdata, handles)
% hObject    handle to exp_vid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile('*.avi', 'Save As...');
if isequal(filename,0) || isequal(pathname,0)
    disp('User selected Cancel');
else
    file = fullfile(pathname,filename);
    fps = str2double(get(handles.fps,'String'));
    movie2avi(handles.Mv,file,'fps',fps);
    disp(['User saved file: ',file]);
end

% --- Executes during object creation, after setting all properties.
function diff_coef_buttons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diff_coef_buttons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function refC_Callback(hObject, eventdata, handles)
% hObject    handle to refC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refC as text
%        str2double(get(hObject,'String')) returns contents of refC as a double


% --- Executes during object creation, after setting all properties.
function refC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_sim.
function run_sim_Callback(hObject, eventdata, handles)
% hObject    handle to run_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

format long;

set(handles.run_sim,'Enable','off');
set(handles.diff_coef_input,'Enable','off');
set(handles.diff_coef_file,'Enable','off');
set(handles.info_text,'String','');
set(handles.run_sim,'Enable','off');
set(handles.exp_fig1,'Enable','off');
set(handles.exp_fig2,'Enable','off');
set(handles.exp_fig3,'Enable','off');
set(handles.exp_gui,'Enable','off');
set(handles.reset,'Enable','on');
set(handles.reset,'Value',get(handles.reset,'Min'));

% Choose figure to export as video
% Best not save frames for all figures or else animation will be too laggy
if isfield(handles,'exp_sel')
    switch handles.exp_sel % only exists if selected object changes
        case 'exp_fig1'
            vid = handles.Conc_3D;
        case 'exp_fig2'
            vid = handles.Conc_Cal;
        case 'exp_fig3'
            vid = handles.Conc_Deriv;
        case 'exp_gui'
            vid = gcf;
    end
else
    vid = handles.Conc_3D; % default action
end

sim_time = str2double(get(handles.tt,'String')); % Simulation time
initC = str2double(get(handles.initC,'String'))*10^str2double(get(handles.pow_initC,'String')); % Initial concentration per unit length at release
diff_coef = str2double(get(handles.diff_coef,'String'))*10^str2double(get(handles.pow_diff_coef,'String')); % diffusion coefficient
refC = str2double(get(handles.refC,'String'))*10^str2double(get(handles.pow_refC,'String')); % ambient glucose concentration
handles.nframes = str2double(get(handles.numframes,'String')); % number of frames

% dimension unit conversion factor to dm for calculation
dim_unit = get(handles.unit,'String');
switch dim_unit
    case 'nm'
        conv_factor = 10^-8;
    case 'um'
        conv_factor = 10^-5;
    case 'mm'
        conv_factor = 10^-2;
    case 'cm'
        conv_factor = 10^-1;
    case 'dm'
        conv_factor = 1;
    case 'm'
        conv_factor = 10;
end

% diffusion dimensions (range)
length_neg = str2double(get(handles.x_neg,'String'));
length_pos = str2double(get(handles.x_pos,'String'));
width_neg = str2double(get(handles.y_neg,'String'));
width_pos = str2double(get(handles.y_pos,'String'));

% Store the electrode coordinates in matrix
% Row number is electrode number, column 1 and 2 are x and y coordinates
elect(1,1) = str2double(get(handles.elect1x,'String'));
elect(1,2) = str2double(get(handles.elect1y,'String'));
elect(2,1) = str2double(get(handles.elect2x,'String'));
elect(2,2) = str2double(get(handles.elect2y,'String'));
elect(3,1) = str2double(get(handles.elect3x,'String'));
elect(3,2) = str2double(get(handles.elect3y,'String'));

% time range
tlist=linspace(0,sim_time,handles.nframes);

% Get the axis interval size (fixed number of intervals)
length_diff = (length_pos - length_neg)/100;
width_diff = (width_pos - width_neg)/100;

% Calculate the corresponding indices in dimensions
elect_ind = round(1+horzcat((elect(1:3,1)-length_neg)/length_diff,(elect(1:3,2)-width_neg)/width_diff));

% Rectangular grid: creates output coordinate arrays x and y so that all
% combinations of x and y can be obtained and used with element arithmetic
[x,y] = meshgrid(length_neg:length_diff:length_pos,width_neg:width_diff:width_pos);

% Convert to dm
x = x*conv_factor;
y = y*conv_factor;
axis_3D = [conv_factor*length_neg conv_factor*length_pos conv_factor*width_neg conv_factor*width_pos refC 1000*initC];

% time-animation generation
% begins here

newplot;
%linkdata on

% pre-allocate for improved speed
electrode = zeros([handles.nframes 3]);
timemat = zeros([handles.nframes 1]);

set(handles.pause_tog,'Enable','on'); % Enables pause/resume toggle button

warning('off','MATLAB:divideByZero'); % hide divide by zero warnings
warning('off','MATLAB:polyfit:PolyNotUnique'); % hide polyfit warnings

% Diffusion in 2D
for n=1:handles.nframes
  
  % z-axis concentration values at each time instant *IMPORTANT*
  Conc = refC + initC/(4*pi*diff_coef*tlist(n))*exp(-(x.^2+y.^2)/(4*diff_coef*tlist(n)));
  
  % Replaces non-finite elements (including NaN) with 0 (must have
  % otherwise cannot display derivatives!)
  i=find(~isfinite(Conc));
  Conc(i)=zeros(size(i));
  
  Concmax = max(max(Conc));
  Concmin = min(min(Conc));

  % Electrodes 1-3 in matrix columns 1-3 respectively
  for nn=1:3
    electrode(n,nn) = Conc(elect_ind(nn,1),elect_ind(nn,2));
  end
  
  timemat(n,:) = tlist(n);
  
  % Due to pre-allocation, take all data up till the loop index n
  l=polyfit(timemat(1:n), electrode(1:n,1),5); % creates polynomial and finds the coefficients that fit the data
  kl=polyder(l); % returns derivative of polynomial
  ll=polyval(kl,timemat(1:n)); % returns evaluation of polynomial at each time instance
  
  p=polyfit(timemat(1:n), electrode(1:n,2),5);
  kp=polyder(p);
  kk=polyval(kp,timemat(1:n));
  
  m=polyfit(timemat(1:n), electrode(1:n,3),5);
  km=polyder(m);
  kmm=polyval(km,timemat(1:n));
  
  % 3D Dispersion Simulation
  surf(handles.Conc_3D,x,y,Conc);
  caxis([Concmin Concmax]);
  colormap(jet),...
  shading(handles.Conc_3D, 'interp');
  axis(handles.Conc_3D,axis_3D);...
      
  % Concentration of Source & Cal-points over time
  plot(handles.Conc_Cal,timemat(1:n), electrode(1:n,1), 'b', timemat(1:n), electrode(1:n,2), 'g', timemat(1:n), electrode(1:n,3), 'r')
  leg=legend(handles.Conc_Cal,'electrode1', 'electrode2','electrode3');
  ylabel(handles.Conc_Cal,'Conc. (mmol/L)');
  xlabel(handles.Conc_Cal,'time');
  grid (handles.Conc_Cal);
  set(leg,'Location', 'NorthWest');
  
  % 1st Derivative of Cal-points
  plot(handles.Conc_Deriv,timemat(1:n),ll,'--b',timemat(1:n),kk,'--g',timemat(1:n),kmm,'--r')
  grid(handles.Conc_Deriv);

  handles.Mv(:,n) = getframe(vid); % new field "Mv" to structure "handles"
  
  % Halts the loop upon pressing down the pause button, and resumes after
  % it is pressed again
  if get(handles.pause_tog,'Value') == get(handles.pause_tog,'Max')
      waitfor(handles.pause_tog,'Value',get(handles.pause_tog,'Min'));
  end
  
  % Execute reset_Callback in loop and return to invoking run_sim_Callback
  if get(handles.reset,'Value') == get(handles.reset,'Max')
      set(handles.reset,'Value',get(handles.reset,'Min'));
      set(handles.reset,'Enable','off');
      return;
  end
  
end

handles.finished = 1; % flag for simulation finished

set(handles.exp_vid,'Enable','on'); % Enables video clip export button
set(handles.pause_tog,'Enable','off');

% Create workspace variables to store data
assignin('base','elect1',electrode(:,1));
assignin('base','elect2',electrode(:,2));
assignin('base','elect3',electrode(:,3));
assignin('base','times',timemat);

guidata(hObject,handles); % save GUI data (change to structure)

function tt_Callback(hObject, eventdata, handles)
% hObject    handle to tt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents  of tt as text
%        str2double(get(hObject,'String')) returns contents of tt as a double


% --- Executes during object creation, after setting all properties.
function tt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numframes_Callback(hObject, eventdata, handles)
% hObject    handle to numframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numframes as text
%        str2double(get(hObject,'String')) returns contents of numframes as a double


% --- Executes during object creation, after setting all properties.
function numframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fps_Callback(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps as text
%        str2double(get(hObject,'String')) returns contents of fps as a double


% --- Executes during object creation, after setting all properties.
function fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in diff_coef_buttons.
function diff_coef_buttons_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in diff_coef_buttons 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag')
    case 'diff_coef_file'
        set(handles.diff_coef,'Enable','off');
        set(handles.pow_diff_coef,'Enable','off');
        [filename, pathname] = uigetfile('*.dat','Select the data file');
    case 'diff_coef_input'
        set(handles.diff_coef,'Enable','on');
        set(handles.diff_coef,'String','1');
        set(handles.pow_diff_coef,'Enable','on');
end


% --- Executes on button press in exp_img.
function exp_img_Callback(hObject, eventdata, handles)
% hObject    handle to exp_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile('*.jpg', 'Save As...');
if isequal(filename,0) || isequal(pathname,0)
    disp('User selected Cancel');
else
    file = fullfile(pathname,filename);
    if isfield(handles,'exp_sel')
        switch handles.exp_sel % only exists if selected object changes
            case 'exp_fig1'
                img = getframe(handles.Conc_3D);
            case 'exp_fig2'
                img = getframe(handles.Conc_Cal);
            case 'exp_fig3'
                img = getframe(handles.Conc_Deriv);
            case 'exp_gui'
                img = getframe(gcf);
        end
    else
        img = getframe(handles.Conc_3D); % default action
    end
    
    if isempty(img.colormap)
        imwrite(img.cdata,file); % RGB images
    else
        imwrite(img.cdata,img.colormap,file); % indexed images
    end
    disp(['User saved file: ',file]);
end


% --- Executes when selected object is changed in exp_buttons.
function exp_buttons_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in exp_buttons 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.exp_sel = get(hObject,'Tag');
guidata(hObject,handles);


% --------------------------------------------------------------------
function exp_buttons_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to exp_buttons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reset

% Clear axes and reset to default state
cla(handles.Conc_3D);
cla(handles.Conc_Cal);
cla(handles.Conc_Deriv);
set(handles.pause_tog,'Value',get(handles.pause_tog,'Min'));
set(handles.pause_tog,'String','Pause');
set(handles.pause_tog,'Enable','off');
set(handles.exp_vid,'Enable','off');
set(handles.run_sim,'Enable','on');
set(handles.exp_fig1,'Enable','on');
set(handles.exp_fig2,'Enable','on');
set(handles.exp_fig3,'Enable','on');
set(handles.exp_gui,'Enable','on');
set(handles.exp_img,'Enable','off');
set(handles.info_text,'String','Choose the figure to export as an video');
set(handles.diff_coef_input,'Enable','on');
set(handles.diff_coef_file,'Enable','on');

if isfield(handles,'finished')
    evalin('base',['clear ' 'elect1 elect2 elect3 times']);
    handles = rmfield(handles,'finished'); % remove field from handles
    guidata(hObject,handles); % need to save to remove field
    set(hObject,'Value',get(hObject,'Min'));
end

function x_neg_Callback(hObject, eventdata, handles)
% hObject    handle to x_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_neg as text
%        str2double(get(hObject,'String')) returns contents of x_neg as a double


% --- Executes during object creation, after setting all properties.
function x_neg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_pos_Callback(hObject, eventdata, handles)
% hObject    handle to x_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_pos as text
%        str2double(get(hObject,'String')) returns contents of x_pos as a double


% --- Executes during object creation, after setting all properties.
function x_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elect1x_Callback(hObject, eventdata, handles)
% hObject    handle to elect1x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elect1x as text
%        str2double(get(hObject,'String')) returns contents of elect1x as a double


% --- Executes during object creation, after setting all properties.
function elect1x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elect1x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elect1y_Callback(hObject, eventdata, handles)
% hObject    handle to elect1y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elect1y as text
%        str2double(get(hObject,'String')) returns contents of elect1y as a double


% --- Executes during object creation, after setting all properties.
function elect1y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elect1y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elect2x_Callback(hObject, eventdata, handles)
% hObject    handle to elect2x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elect2x as text
%        str2double(get(hObject,'String')) returns contents of elect2x as a double


% --- Executes during object creation, after setting all properties.
function elect2x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elect2x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elect2y_Callback(hObject, eventdata, handles)
% hObject    handle to elect2y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elect2y as text
%        str2double(get(hObject,'String')) returns contents of elect2y as a double


% --- Executes during object creation, after setting all properties.
function elect2y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elect2y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elect3x_Callback(hObject, eventdata, handles)
% hObject    handle to elect3x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elect3x as text
%        str2double(get(hObject,'String')) returns contents of elect3x as a double


% --- Executes during object creation, after setting all properties.
function elect3x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elect3x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function elect3y_Callback(hObject, eventdata, handles)
% hObject    handle to elect3y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of elect3y as text
%        str2double(get(hObject,'String')) returns contents of elect3y as a double


% --- Executes during object creation, after setting all properties.
function elect3y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to elect3y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_neg_Callback(hObject, eventdata, handles)
% hObject    handle to y_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_neg as text
%        str2double(get(hObject,'String')) returns contents of y_neg as a double


% --- Executes during object creation, after setting all properties.
function y_neg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_pos_Callback(hObject, eventdata, handles)
% hObject    handle to y_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_pos as text
%        str2double(get(hObject,'String')) returns contents of y_pos as a double


% --- Executes during object creation, after setting all properties.
function y_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pow_initC_Callback(hObject, eventdata, handles)
% hObject    handle to pow_initC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pow_initC as text
%        str2double(get(hObject,'String')) returns contents of pow_initC as a double


% --- Executes during object creation, after setting all properties.
function pow_initC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pow_initC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pow_refC_Callback(hObject, eventdata, handles)
% hObject    handle to pow_refC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pow_refC as text
%        str2double(get(hObject,'String')) returns contents of pow_refC as a double


% --- Executes during object creation, after setting all properties.
function pow_refC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pow_refC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pow_diff_coef_Callback(hObject, eventdata, handles)
% hObject    handle to pow_diff_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pow_diff_coef as text
%        str2double(get(hObject,'String')) returns contents of pow_diff_coef as a double


% --- Executes during object creation, after setting all properties.
function pow_diff_coef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pow_diff_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit_Callback(hObject, eventdata, handles)
% hObject    handle to unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit as text
%        str2double(get(hObject,'String')) returns contents of unit as a double


% --- Executes during object creation, after setting all properties.
function unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
