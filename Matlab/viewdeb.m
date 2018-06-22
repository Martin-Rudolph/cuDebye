function viewdeb(File,Limits)
% Plot structure from deb file. As this function is very slow for large
% clusters, limits for plotting are set commonly from -10 to 10.
%
% You can specify limits:
% viewdeb([-50, 50]) or viewdeb([-50, 50; -20, 50; -10, 10])
%
% You can specify also a File:
% viewdeb(Filename)
% 
% Or both: viewdeb(Filename,Limits)


    if nargin==0
        [FileName, PathName] = uigetfile('*.deb','Select a File for Plotting');
        File = [PathName,FileName];
        Limits = [-10, 10];
    end
    if nargin==1 && ischar(File)==0
        Limits = File;
        [FileName, PathName] = uigetfile('*.deb','Select a File for Plotting');
        File = [PathName,FileName];
    elseif nargin==1 && ischar(File)==1
        Limits = [-10, 10];
    end

    fid = fopen(File,'r');
    for i=1:4
        Line = fgetl(fid);
    end
    nSG = str2double(Line);
    for i=1:nSG
        Line = fgetl(fid);
        Line = strsplit(Line);
        nAt(i) = str2double(Line{1});
        AtomType{i} = Line{2};
    end
    fclose(fid);
    Data = dlmread(File,'',4+nSG,0);

    Data = [ones(size(Data,1),1),Data];
    Start = [0,cumsum(nAt)]+1;
    ElementNumber = findelement(AtomType);
    for i=1:nSG
        Data(Start(i):Start(i+1)-1,1) = Data(Start(i):Start(i+1)-1,1)*ElementNumber(i);
    end
    
    if isempty(Limits)==0
        if size(Limits,1)==1
            Limits = repmat(Limits,3,1);
        end
        for i=1:3
            Data = Data(Data(:,i+1)>Limits(i,1) & Data(:,i+1)<Limits(i,2),:);
        end
    end

    Fig = figure('Name',FileName);
    showstruct(Data,Fig)

end

function  showstruct( Data, Fig )
% Plot structure
    if nargin==1
        Fig = figure(1);clf;
    end
    MarkerSize = 5;
    Cell = false;
    
    ElementNumber = Data(:,1);
    x = Data(:,2);
    y = Data(:,3);
    z = Data(:,4);
    
    uElements  = unique(ElementNumber);
    [~,uColor] = getradii(uElements);
    [~,uElements]  = findelement(uElements);
    
    [Radii, Color] = getradii(ElementNumber);
    A              = pi/4*(2*Radii*MarkerSize).^2;

    Fig; hold on
    scatter3(x,y,z,A,Color, 'filled', 'MarkerEdgeColor','k')
    xlabel('x [Å]')
    ylabel('y [Å]')
    zlabel('z [Å]')
    for i=1:length(uElements)
        text(1.1,1-(i-1)*0.1, uElements{i}, 'Units','normalized', 'Color',uColor(i,:), 'Fontsize',24, 'FontWeight','bold');
    end
    
    set(gca, 'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio',[1 1 1],...
        'Xlim', [min(x),max(x)], 'Ylim', [min(y),max(y)], 'Zlim', [min(z),max(z)])
    hold off

end

function [ElementNumber, ElementSymbol] = findelement(Element)
    TP = false;
    if size(Element,2)>size(Element,1)
        Element = Element';
        TP = true;
    end
    

    FunctionPath = mfilename('fullpath');
    Index = strfind(FunctionPath,'\');
    FunctionPath = FunctionPath(1:Index(end));
    fid   = fopen([FunctionPath,'Tables\elements.txt'],'r');
    Data  = textscan(fid,'%f%s%s');
    fclose(fid);
    
    
    Number = Data{1};
    Symbol = Data{2};
    Name   = Data{3};
    
    if iscell(Element)
        Element = regexprep(Element,'-\w*','');
        Element = regexprep(Element,'+\w*','');
        Element = regexprep(Element,'(\w*','');
        Element = regexprep(Element,':\w*','');
        Element = regexprep(Element,'\d\w*','');
    end
    
    UniqueType = unique(Element);
    
    if iscell(Element)
        ElementNumber = zeros(size(Element));
        ElementSymbol = Element;
        for i=1:length(UniqueType)
            Index  = cellfun(@(x) strcmpi(x,UniqueType{i}),Symbol);
            Num = Number(Index); %find(Index==1);
            Index = cellfun(@(x) strcmpi(x,UniqueType{i}),Element);
            if isempty(Num)==0
                ElementNumber(Index) = Num*ones(sum(Index),1);
            end
        end    
    elseif isnumeric(Element)
        ElementSymbol = cell(size(Element));
        ElementNumber = Element;
        for i=1:length(UniqueType)
            Symb = Symbol(Number(Number==UniqueType(i)));
            ElementSymbol(Element==UniqueType(i)) = Symb;
        end
    end
    
    if TP==true
        ElementNumber = ElementNumber';
        ElementSymbol = ElementSymbol';
    end
    
end

function [ElementRadii, ElementColor] = getradii( Element )
    TP = false;
    if size(Element,2)>size(Element,1)
        Element = Element';
        TP = true;
    end

    FunctionPath = mfilename('fullpath');
    Index = strfind(FunctionPath,'\');
    FunctionPath = FunctionPath(1:Index(end));
    fid   = fopen([FunctionPath,'Tables\atomicradii.txt'],'r');
    Data  = textscan(fid,'%s%f%f%f%f');
    fclose(fid);
    
    Symbol = Data{:,1};
    Radii  = Data{:,2};
    RGB    = cell2mat(Data(:,3:5)); 
    
    if isnumeric(Element)
        [~,Element] = findelement(Element);
    end
    
    ElementRadii = zeros(size(Element));
    ElementColor = zeros(size(Element,1),3);
    
    UniqueType = unique(Element);
    for i=1:length(UniqueType)
        Index  = cellfun(@(x) strcmpi(x,UniqueType{i}),Symbol);
        R = Radii(Index);
        Col = RGB(Index,:);
        Index = cellfun(@(x) strcmpi(x,UniqueType{i}),Element);
        Ones = ones(sum(Index),1);
        ElementRadii(Index) = R*Ones;
        ElementColor(Index,1) = Col(1)*Ones;
        ElementColor(Index,2) = Col(2)*Ones;
        ElementColor(Index,3) = Col(3)*Ones;
    end
    
end


