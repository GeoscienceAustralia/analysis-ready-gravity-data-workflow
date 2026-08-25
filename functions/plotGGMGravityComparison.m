function plotGGMGravityComparison( ...
    Long, Lat, GravAnom, GGMinterpolated, outputFile)
% plotWAGravityComparison
% Plots WA gravity anomaly and (WA – EGM2008) difference,
% computes the mean difference, and saves the figure.
%
% Inputs:
%   Long              - Longitude vector
%   Lat               - Latitude vector
%   GravAnom          - WA gravity anomalies (mGal)
%   GGMinterpolated   - Interpolated EGM2008 anomalies (mGal)
%   outputFile        - Output image filename (e.g. 'outputs/plots/EGM2008WA.png')
%
% Output:
%   meanDiff              - Mean difference (WA – EGM2008) in mGal

    figure
    hold on
    scatter(Long, Lat, 1, GravAnom)
    colormap(jet)
    cb1 = colorbar;
    title(cb1,'mGal','FontSize',10)
    title('Gravity anomaly')

    figure
    hold on
    diffVals = GravAnom - GGMinterpolated;
    scatter(Long, Lat, 1, diffVals)
    colormap(jet)
    cb2 = colorbar;
    title(cb2,'mGal','FontSize',10)
    title('EGM2008 - gravity anomaly')

    % ---------- Compute mean ----------
    meanDiff = mean(diffVals,'omitnan');

    % ---------- Add mean as figure text ----------
    annotation('textbox',[0.15 0.02 0.7 0.05], ...
        'String', sprintf('Mean(EGM2008 – gravity anomaly): %.4f mGal', meanDiff), ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'FontSize',10, ...
        'FontWeight','bold');

    % ---------- Save ----------
    if nargin >= 5 && ~isempty(outputFile)
        saveas(gcf, outputFile);
    end
end