function plotdebi(Normalized)
% Plots multiple debI-files and the resulting average.
% With the additional input arguments 'int' and 'max' the intensities can
% be normed to integral or maximal intentsity.
    
    [FileName, PathName] = uigetfile('*.debI','MultiSelect','on');
    if nargin==0
        Normalized = 'none';
    end

    if ischar(FileName)
        FileName = {FileName};
    end
    Xmin = 0;
    Xmax = 180;
    for i=1:length(FileName)
        DataTemp = dlmread([PathName,FileName{i}],'\t',0,0);
        Data(i,1) = {DataTemp};
        X = Data{i,1}(:,1);
        Y = Data{i,1}(:,2);
        if min(X)>Xmin
            Xmin = min(X);
        end
        if max(X)<Xmax
            Xmax = max(X);
        end
    end
    for i=1:length(FileName)
        DataTemp = dlmread([PathName,FileName{i}],'\t',0,0);
        Data(i,1) = {DataTemp};
        X = Data{i,1}(:,1);
        Y = Data{i,1}(:,2);
        Idx = find(X>=Xmin & X<=Xmax);
        switch lower(Normalized)
            case {'int','integral'}
                Data(i,2) = {[X,Y./trapz(X(Idx),Y(Idx))]};
                LabelY = 'Normalized Intensity (I/Iint)';
            case {'max','maximum','maximal'}
                Data(i,2) = {[X,Y./max(Y(Idx))]};
                LabelY = 'Normalized Intensity (I/Imax)';
            otherwise
                Data(i,2) = {[X,Y]};
                LabelY = 'Intensity [e.u.]';
        end
    end

    aTT = (Xmin:0.02:Xmax)';
    aI  = zeros(size(aTT));
    figure(1);clf
    set(gca,'xLim',[Xmin,Xmax])
    hold on
    for i=1:length(FileName)
        plot(Data{i,2}(:,1),Data{i,2}(:,2))
        aI = aI + spline(Data{i,1}(:,1),Data{i,1}(:,2),aTT);
    end
    switch lower(Normalized)
        case {'int','integral'}
            aI = aI/trapz(aTT,aI);
        case {'max','maximum','maximal'}
            aI = aI/max(aI);
        otherwise
            aI = aI/length(FileName);
    end
    plot(aTT,aI,'k', 'LineWidth',2)
    legend([FileName,'Average'], 'Interpreter','none')
    xlabel('Diffraction angle [°]')
    ylabel(LabelY)
    hold off
end