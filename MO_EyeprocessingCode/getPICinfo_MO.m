function [picINFOtab] = getPICinfo_MO(inTable)


picJPG = inTable.ClipName;
picPARTS = split(picJPG,{'_','.'});
picINFOtab = inTable;
picINFOtab.BorderType = picPARTS(:,1);
picINFOtab.CutType = picPARTS(:,2);
picINFOtab.ClipNUM = picPARTS(:,3);

% tmpCombine = cell(length(catSubn),1);
% for pi = 1:length(catSubn)
%     tmpCombine{pi} = [num2str(catSubn(pi)),'.',num2str(picNUM(pi))];
% end
% picINFOtab.CatPICid = tmpCombine;

end