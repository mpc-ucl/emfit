
function [D,Z] = getData(dataDir,selectors)

    dr   = dir(dataDir);
    count = 1;
    for i=1:length(dr)
        doFile = true;
        for j=1:length(selectors)
            if isempty(strfind(dr(i).name, selectors{j}))
                doFile = false;
            end
        end
       if doFile
           load([dataDir, dr(i).name]);
           results.id = dr(i).name(1:10);
           D(count) = results;
           Z = params;
           count = count + 1;
       end
    end
    
    [D,Z] = precomputeParams(D,Z);

end