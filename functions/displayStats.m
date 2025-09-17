
function displayStats(arr)
    % displayStats - Displays basic statistics for a numeric array
    % Usage: displayStats(arr)
    % Example: displayStats([1, 2, 3, 4, 5])

    if ~isnumeric(arr)
        error('Input must be a numeric array.');
    end

    fprintf('Statistics for the input array:\n');
    fprintf('--------------------------------\n');
    fprintf('Minimum: %.2f\n', round(min(arr)));
    fprintf('Maximum: %.2f\n', round(max(arr)));
    fprintf('Mean: %.2f\n',    round(mean(arr, 'omitnan')));
    fprintf('Median: %.2f\n',  round(median(arr, 'omitnan')));
    fprintf('Standard Deviation: %.2f\n', round(std(arr, 'omitnan')));
    %fprintf('Variance: %.2f\n', var(arr, 'omitnan'));
    %fprintf('Sum: %.2f\n', sum(arr, 'omitnan'));
    fprintf('Length: %d\n', length(arr));
end


