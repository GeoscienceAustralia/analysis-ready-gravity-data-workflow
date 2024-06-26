function varargout = custom_grpstats(data, groupingVar, functionsCell)
    % CUSTOM_GRPSTATS Computes grouped statistics for specified functions.
    %
    %   VARARGOUT = CUSTOM_GRPSTATS(DATA, GROUPINGVAR, FUNCTIONCELL) computes
    %   grouped statistics for the data in DATA based on the grouping
    %   specified by GROUPINGVAR. FUNCTIONCELL is a cell array of function
    %   handles specifying the statistics to be calculated for each group.
    %
    %   VARARGOUT contains separate output variables for each function in
    %   FUNCTIONCELL, where each variable corresponds to a specific
    %   statistic for each group.
    %
    %   Example:
    %   data = [1; 2; 3; 4; 5; 6];
    %   groupingVar = [1; 1; 2; 2; 3; 3];
    %   functionsCell = {@sum, @mean, @std};  % Example: sum, mean, and standard deviation
    %
    %   [sumStats, meanStats, stdStats] = custom_grpstats(data, groupingVar, functionsCell);
    %
    % Written by ChatGPT3.5
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    uniqueGroups = unique(groupingVar);

    numGroups = numel(uniqueGroups);
    numStats = numel(functionsCell);

    varargout = cell(1, numStats);

    for j = 1:numStats
        varargout{j} = zeros(numGroups, 1);
    end

    for i = 1:numGroups
        groupIdx = (groupingVar == uniqueGroups(i));
        groupData = data(groupIdx);

        for j = 1:numStats
            currentFunction = functionsCell{j};
            varargout{j}(i) = feval(currentFunction, groupData);
        end
    end
end


