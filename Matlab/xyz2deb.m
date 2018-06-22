function xyz2deb(File,Prop)
% Keep the function empty for GUI -> xyz2deb
% Pass a string for file directory and name to use function without GUI, 
% e.g. in a matlab script.
% 
% As second argument you can also pass a structure variable to specify the
% header parameters in the deb-file. For more information pass a random 
% parameter first and check the output.
%
% Alternatively you can edit the standard properties within the function.

%% Here you can enter your default Cuda device and diffraction properties
    if nargin<2 
        Device = 0;
        BlockSize = 16;
        GridSize = 2048;

        Wavelength = 1.5418;
        TwoThetaMin = 10;
        TwoThetaStep = 0.02;
        TwoThetaMax = 140;

        Rstep = 0.0001;
        
        I = true;
        R = false;
    %% Do not edit after this line
    elseif nargin==2
        try
            Device    = Prop.Cuda.Device;
            BlockSize = Prop.Cuda.BlockDimXY;
            GridSize  = Prop.Cuda.GridDimXY;

            Wavelength = Prop.Diff.Wavelength;
            TwoThetaMin = Prop.Diff.TwoThetaMin;
            TwoThetaStep = Prop.Diff.TwoThetaStep;
            TwoThetaMax = Prop.Diff.TwoThetaMax;
            
            try
                I = Prop.Out.I;
                R = Prop.Out.R;
            catch
                I = true;
                R = false;
            end

            Rstep = Prop.Hist.Rstep;
        catch
            disp('The second argument has to be a single structure variable containing the following information:')
            disp('Prop.Cuda.Device')
            disp('Prop.Cuda.BlockDimXY')
            disp('Prop.Cuda.GridDimXY')
            disp('Prop.Diff.Wavelength')
            disp('Prop.Diff.TwoThetaMin')
            disp('Prop.Diff.TwoThetaStep')
            disp('Prop.Diff.TwoThetaMax')
            disp('Prop.Out.I -> true or false')
            disp('Prop.Out.R -> true or false')
        end
    else
        disp('ERROR: Too many input arguments!')
    end
    if nargin==0
        Handles.GUI = figure('Units','normalized');
        Handles.Text.Device = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.Device.Position = [0.05, 0.9, 0.2, 0.05];
        Handles.Text.Device.String   = 'Cuda Device Index: ';
        Handles.Edit.Device = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.Device.Position = [0.25, 0.9, 0.2, 0.05];
        Handles.Edit.Device.String   = num2str(Device);

        Handles.Text.BlockSize = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.BlockSize.Position = [0.05, 0.8, 0.2, 0.05];
        Handles.Text.BlockSize.String   = 'Block Size: ';
        Handles.Edit.BlockSize = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.BlockSize.Position = [0.25, 0.8, 0.2, 0.05];
        Handles.Edit.BlockSize.String   = num2str(BlockSize);

        Handles.Text.GridSize = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.GridSize.Position = [0.05, 0.7, 0.2, 0.05];
        Handles.Text.GridSize.String   = 'Grid Size: ';
        Handles.Edit.GridSize = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.GridSize.Position = [0.25, 0.7, 0.2, 0.05];
        Handles.Edit.GridSize.String   = num2str(GridSize);

        Handles.Text.Rstep = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.Rstep.Position = [0.05, 0.6, 0.2, 0.05];
        Handles.Text.Rstep.String   = 'Rstep [A]: ';
        Handles.Edit.Rstep = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.Rstep.Position = [0.25, 0.6, 0.2, 0.05];
        Handles.Edit.Rstep.String   = num2str(Rstep);

        Handles.Text.Wavelength = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.Wavelength.Position = [0.5, 0.9, 0.2, 0.05];
        Handles.Text.Wavelength.String   = 'Wavelength [A]: ';
        Handles.Edit.Wavelength = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.Wavelength.Position = [0.7, 0.9, 0.2, 0.05];
        Handles.Edit.Wavelength.String   = num2str(Wavelength);

        Handles.Text.TwoThetaMin = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.TwoThetaMin.Position = [0.5, 0.8, 0.2, 0.05];
        Handles.Text.TwoThetaMin.String   = 'TwoThetaMin [°]: ';
        Handles.Edit.TwoThetaMin = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.TwoThetaMin.Position = [0.7, 0.8, 0.2, 0.05];
        Handles.Edit.TwoThetaMin.String   = num2str(TwoThetaMin);

        Handles.Text.TwoThetaMax = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.TwoThetaMax.Position = [0.5, 0.7, 0.2, 0.05];
        Handles.Text.TwoThetaMax.String   = 'TwoThetaMax [°]: ';
        Handles.Edit.TwoThetaMax = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.TwoThetaMax.Position = [0.7, 0.7, 0.2, 0.05];
        Handles.Edit.TwoThetaMax.String   = num2str(TwoThetaMax);

        Handles.Text.TwoThetaStep = uicontrol('Parent',Handles.GUI, 'Style','text', 'Units','normalized');
        Handles.Text.TwoThetaStep.Position = [0.5, 0.6, 0.2, 0.05];
        Handles.Text.TwoThetaStep.String   = 'TwoThetaStep [°]: ';
        Handles.Edit.TwoThetaStep = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.TwoThetaStep.Position = [0.7, 0.6, 0.2, 0.05];
        Handles.Edit.TwoThetaStep.String   = num2str(TwoThetaStep);
        
        Handles.RadioButton.I = uicontrol('Parent',Handles.GUI, 'Style','radiobutton', 'Units','normalized');
        Handles.RadioButton.I.Position = [0.7, 0.5, 0.2, 0.05];
        Handles.RadioButton.I.String   = 'Intensity File';
        Handles.RadioButton.I.Value    = I;
        Handles.RadioButton.R = uicontrol('Parent',Handles.GUI, 'Style','radiobutton', 'Units','normalized');
        Handles.RadioButton.R.Position = [0.25, 0.5, 0.2, 0.05];
        Handles.RadioButton.R.String   = 'Distance Files';
        Handles.RadioButton.R.Value    = R;
        

        Handles.Edit.File = uicontrol('Parent',Handles.GUI, 'Style','Edit', 'Units','normalized');
        Handles.Edit.File.Position = [0.05, 0.4, 0.9, 0.05];
        Handles.Edit.File.String   = '';
        Handles.PushButton.File = uicontrol('Parent',Handles.GUI, 'Style','Pushbutton', 'Units','normalized');
        Handles.PushButton.File.Position = [0.05, 0.275, 0.9, 0.1];
        Handles.PushButton.File.String   = 'Browse';
        Handles.PushButton.File.Callback = {@GetFile,Handles};

        Handles.PushButton.Ok = uicontrol('Parent',Handles.GUI, 'Style','pushbutton', 'Units','normalized');
        Handles.PushButton.Ok.Position = [0.4, 0.05, 0.2, 0.075];
        Handles.PushButton.Ok.String   = 'Ok';
        Handles.PushButton.Ok.Callback = {@xyztodeb,Handles};
    end
    if nargin>0
       xyztodeb(0,0,0,File,Device,BlockSize,GridSize,Wavelength,TwoThetaMin,TwoThetaMax,TwoThetaStep,Rstep,I,R)
    end
end

function GetFile(~,~,Handles)
    [FileName, PathName] = uigetfile('*.xyz','Pick a xyz File!');
    Handles.Edit.File.String = [PathName,FileName];
end

function xyztodeb(~,~,Handles,File,Device,BlockSize,GridSize,Wavelength,TwoThetaMin,TwoThetaMax,TwoThetaStep,Rstep,I,R)
    if nargin==3
        File = Handles.Edit.File.String;
    end
    fid = fopen(File);
    Line = fgetl(fid);
    Line = strsplit(Line);
    Columns = 0;
    for i=1:length(Line)
        if isempty(Line{i})==0
            Columns = Columns+1;
        end
    end
    fclose(fid);
    ImportString = ['%s',repmat('%f',1,Columns-1)];

    fid = fopen(File,'r');
    Data = textscan(fid,ImportString);
    fclose(fid);
    Element = Data{1};
    UnElement = unique(Element);
    AtomTyp = zeros(size(Element));
    for i=1:length(UnElement)
        BoolIndex = cellfun(@(String) strcmp(String,UnElement{i}), Data{1});
        AtomTyp = AtomTyp+BoolIndex*i;
    end
    if Columns==4
        Data = [AtomTyp, Data{2:end}, ones(length(AtomTyp),1), zeros(length(AtomTyp),1)];
    elseif Columns==5
        Data = [AtomTyp, Data{2:end}, zeros(length(AtomTyp),1)];
    elseif Columns==6
        Data = [AtomTyp, Data{2:end}];
    end
    assignin('base','Data',Data)
    
    Data  = sortrows(Data,[1,5,6]);
    Index = Data(2:end,1)-Data(1:end-1,1)+Data(2:end,5)-Data(1:end-1,5)+Data(2:end,6)-Data(1:end-1,6);
    Index = find(Index~=0);
    Index = [0;Index;size(Data,1)];
    Start = Index(1:end-1)+1;
    End   = Index(2:end);
    
    D = End-Start+1;
    
    for i=2:4
        Max = max(Data(:,i));
        Min = min(Data(:,i));
        Data(:,i) = Data(:,i)-Min-(Max-Min)/2;
    end
    [~,~,r] = cart2sph(Data(:,2),Data(:,3),Data(:,4));
    Rmax = 2*ceil(max(r)*1.025);
    
    if nargin==3
        Device = str2double(Handles.Edit.Device.String);
        BlockSize = str2double(Handles.Edit.BlockSize.String);
        GridSize = str2double(Handles.Edit.GridSize.String);
        Rstep = str2double(Handles.Edit.Rstep.String);

        Wavelength   = str2double(Handles.Edit.Wavelength.String);
        TwoThetaMin  = str2double(Handles.Edit.TwoThetaMin.String);
        TwoThetaMax  = str2double(Handles.Edit.TwoThetaMax.String);
        TwoThetaStep = str2double(Handles.Edit.TwoThetaStep.String);

        I = Handles.RadioButton.I.Value;
        R = Handles.RadioButton.R.Value;
    end
    
    if I==0 && R==0
        Output = 0;
    end
    if I==1 && R==0
        Output = 0;
    end
    if I==0 && R==1
        Output = 1;
    end
    if I==1 && R==1
        Output = 2;
    end
    
    tic    
    File = [File(1:end-3),'deb'];
    fid = fopen(File,'w');
    fprintf(fid, '%0.0f %0.0f %0.0f',Device,BlockSize, GridSize);
    fprintf(fid, '\n%0.6f %0.2f %0.2f %0.2f', Wavelength, TwoThetaMin, TwoThetaMax, TwoThetaStep);
    fprintf(fid, '\n%0.0f %0.4f %0.0f', Rmax, Rstep, Output);
    fprintf(fid, '\n%0.0f',length(D));
    for i=1:length(D)
        nAtoms = D(i);
        Element = UnElement{Data(Start(i),1)};
        Occupancy = Data(Start(i),5);
        Beq = Data(Start(i),6);
        fprintf(fid, '\n%0.0f',nAtoms);
        fprintf(fid, ' %s',Element);
        fprintf(fid, ' %0.4f',Occupancy);
        fprintf(fid, ' %0.4f',Beq);
        [~,a,b,c] = atomicscatter(Element);
        for j=1:5
            fprintf(fid, ' %0.6f',a(j));
        end
        for j=1:5
            fprintf(fid, ' %0.6f',b(j));
        end
        fprintf(fid, ' %0.6f',c);
    end
    fprintf(fid,'\n');
    fclose(fid);
    disp('Writing Atom Positions...')
    dlmwrite(File,Data(:,2:4),'delimiter','\t','-append','precision',7);
    Time = toc;
    
    disp([num2str(Time),'s for creating a deb-file'])
end

function [f,ax,bx,cx] = atomicscatter(Element,Wavelength,TwoTheta)
%ATOMICSCATTER calculates atomic scattering factor and returns the 11 parameter approach 
%  [f,ax,bx,cx] = atomicscatter(Element,Wavelength,TwoTheta)
%  [f,ax,bx,cx] = atomicscatter(Element,s)
%  [~,ax,bx,cx] = atomicscatter(Element)
s = -1;
TransCheck = false;
if nargin==3
    s  = 2*sind(TwoTheta/2)/Wavelength;
elseif nargin==2
    s = Wavelength;
end
if size(s,1)==1
    s = s';
    TransCheck = true;
end

    FunctionPath = mfilename('fullpath');
    Index = strfind(FunctionPath,'\');
    FunctionPath = FunctionPath(1:Index(end));
    fid   = fopen([FunctionPath,'Tables\atomic_scattering_factors_RIETAN.txt'],'r');
    Line = '#';
    while Line(1)=='#'
        Line  = strtrim(fgetl(fid));
    end
    Data = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f');
    fclose(fid);
    Index = cellfun(@(String) strcmpi(String,Element), Data{1});
    Index = find(Index==1);
    Data = [Data{2:end}];
    ax = Data(Index,1:5);
    bx = Data(Index,7:11);
    cx = Data(Index,6);
    
    if nargin==1
        f = [];
        return;
    end

    f = sum((ones(size(s))*ax).*exp(-0.25*power(s,2)*bx),2)+cx;
if TransCheck==true
    f = f';
end
end