function [children]=findchildren(parents,annotated,pathId)
children=[];
for i = 1:numel(pathId)
    if any(ismember(pathId{i},parents))
        children=[children,annotated.id(i)];
    end
end
end
