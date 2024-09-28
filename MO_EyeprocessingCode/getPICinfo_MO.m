function [picINFOtab] = getPICinfo_MO(inTable,mo_sessionID)


switch mo_sessionID
    case 'en'
        picJPG = inTable.ClipName;
        picPARTS = split(picJPG,{'_','.'});
        picINFOtab = inTable;
        picINFOtab.BorderType = picPARTS(:,1);
        picINFOtab.CutType = picPARTS(:,2);
        picINFOtab.ClipNUM = picPARTS(:,3);

    case 'sc'
        borderTYPE = cell(height(inTable),1);
        cutTYPE = cell(height(inTable),1);
        clipTYPE = cell(height(inTable),1);
        picJPG = inTable.ClipName;
        picINFOtab = inTable;
        for ipi = 1:length(picJPG)
            tmpSPlit = split(picJPG{ipi},{'_','.'});
            borderTYPE{ipi,1} = tmpSPlit{1};
            cutTYPE{ipi,1} = tmpSPlit{2};
            clipTYPE{ipi,1} = tmpSPlit{3};
        end
        picINFOtab.BorderType = borderTYPE;
        picINFOtab.CutType = cutTYPE;
        picINFOtab.ClipNUM = clipTYPE;

    case 'ti'
        picJPG = inTable.ClipName;
        picPARTS = split(picJPG,{'_','.'});
        picINFOtab = inTable;
        picINFOtab.BorderType = picPARTS(:,1);
        picINFOtab.CutType = picPARTS(:,2);
        picINFOtab.ClipNUM = picPARTS(:,3);

end

% tmpCombine = cell(length(catSubn),1);
% for pi = 1:length(catSubn)
%     tmpCombine{pi} = [num2str(catSubn(pi)),'.',num2str(picNUM(pi))];
% end
% picINFOtab.CatPICid = tmpCombine;

end