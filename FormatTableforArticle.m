clc
clear
close all

Task = 'Li';
uiopen(['J:\Piano_Fatigue\Data_Exported\TableEvolution_' Task '.xlsx'],1)

TableforArticle = TableEvolutionLi;
TableforArticle = {};
seq = 2:2:15;
for iR = 1:size(TableEvolutionLi,1)
    for iC = seq
        if table2array(TableEvolutionLi(iR,iC)) == 0
            TableforArticle{iR,iC} = {''};
        elseif table2array(TableEvolutionLi(iR,iC)) < 0 & table2array(TableEvolutionLi(iR,iC+1)) < 0.05
            TableforArticle{iR,iC} = '-';
            TableforArticle{iR,iC+1} = -1;
        elseif table2array(TableEvolutionLi(iR,iC)) > 0 & table2array(TableEvolutionLi(iR,iC+1)) < 0.05
            TableforArticle{iR,iC} = '+';
            TableforArticle{iR,iC+1} = -1;
        end
    end
end



