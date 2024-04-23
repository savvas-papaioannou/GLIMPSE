function clfall
    FigList = findall(groot, 'Type', 'figure');
    for iFig = 1:numel(FigList)
        try
            clf(FigList(iFig));
        catch
            % Nothing to do
        end
    end
end