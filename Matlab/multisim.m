function multisim
% Run function for simulation of all structures within a user specified
% folder. For a recalculation you have to delete the calculated debI file.
% If the folder contains a xyz-file but not the corresponding deb-file,
% this file will be created automatically. Be carefull the standard
% parameters from the function xyz2deb will be used.
% In contrast to the other scripts this will work only, when the
% branch folder structure stays intact.
    
    cuDebyePath = pwd;
    Idx = strfind(cuDebyePath,'Matlab');
    cuDebyePath = cuDebyePath(1:Idx(end)-1);

    PathName = uigetdir(cuDebyePath);

    FileNameXYZ = getfiles(PathName,'.xyz');
    FileNameDEB = getfiles(PathName,'.deb');

    FileName = {};
    k = 0;
    for i=1:length(FileNameXYZ)
        FileDEB = FileNameXYZ{i}; 
        FileDEB = [FileDEB(1:end-3),'deb'];
        BoolIndex = cellfun(@(String) strcmp(String,FileDEB), FileNameDEB);
        if sum(BoolIndex)==0
            k = k+1;
            FileName{k} = FileNameXYZ{i};
        end
    end

    for i=1:length(FileName)
        fprintf('Convert File %0.0f of %0.0f -- %s\n',i,length(FileName),FileName{i});
        File = [PathName,'\',FileName{i}];
        xyz2deb(File);
    end

    %%
    FileNameDEB = getfiles(PathName,'.deb');
    FileNameDEBI = getfiles(PathName,'.debI');

    FileName = {};
    k = 0;
    for i=1:length(FileNameDEB)
        FileDEBI = FileNameDEB{i}; 
        FileDEBI = [FileDEBI,'I'];
        BoolIndex = cellfun(@(String) strcmp(String,FileDEBI), FileNameDEBI);
        if sum(BoolIndex)==0
            k = k+1;
            FileName{k} = FileNameDEB{i};
        end
    end
    
    for i=1:length(FileName)
        fprintf('Calculate File %0.0f of %0.0f -- %s\n',i,length(FileName),FileName{i});
        File = [PathName,'\',FileName{i}];
        [Status,cmdOut] = system([cuDebyePath, 'cuDebye.exe ',File],'-echo');
        figure(1);clf;
            hold on
            title(['Actual Calculation: ',regexprep(FileName{i},'_','\\_')]);
            xlabel('Diffraction Angle [°]')
            ylabel('Intensity [e.U.]')
            TempData = dlmread([File,'I'],'\t',0,0);
            plot(TempData(:,1),TempData(:,2))
            hold off
    end

    FileName = getfiles(PathName,'debI');
    figure(1);clf;
        xlabel('Diffraction Angle [°]')
        ylabel('Normalized Intensity I/Imax')
    for i=1:length(FileName)
        File = [PathName,'\',FileName{i}];
        figure(1)
            hold on
            TempData = dlmread(File,'\t',0,0);
            plot(TempData(:,1),TempData(:,2)/max(TempData(:,2)))
            hold off
    end
    legend(FileName,'Interpreter','none')
end

function [ FileName ] = getfiles( ParentPath, varargin )
%getfiles: Returns all filenames of a specific directory and typ.
%   [FileName] = getfiles(ParentPath) returns all files from a Path
%   [FileName] = getfiles(ParentPath,'.dat',...,'.txt') returns all files of a specific typ

    AllFiles = dir(ParentPath);
    IsFile = cell2mat({AllFiles.isdir})==0;
    FileName = {AllFiles(IsFile).name};
    
    if sum(IsFile) == 0
        FileName = [];
        return
    end
    
    if nargin>1
        DataTyp = cell(nargin-1);
        for i = 1:nargin-1;
            DataTyp{i} = varargin{i};
        end
        DeleteFiles = true(size(FileName));
        for i = 1:length(DataTyp)
            typ = DataTyp{i};
            len = length(typ);
            for j = 1:length(FileName)
                file = FileName{j}; file = file(end-len+1:end);
                if strcmp(file,typ)==1
                    DeleteFiles(j)=0;
                end
            end
        end
        FileName(DeleteFiles) = [];
    end
    
end