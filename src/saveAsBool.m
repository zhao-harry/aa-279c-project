function saveAsBool(fig, fileName, savePlots)
    if savePlots == true
        saveas(fig, fileName)
    end
end