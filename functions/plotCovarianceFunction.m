function plotCovarianceFunction(gravityData, residuals, outputParameters, gravityType, covarianceInfo, covarianceParameters, fittedCovariance,block_counter)

    MINLONG=min(gravityData(:,1));
    MAXLONG=max(gravityData(:,1));
    MINLAT=min(gravityData(:,2));
    MAXLAT=max(gravityData(:,2));
    % for writing covariance parameters on the plots
    str = {['parameter A ',num2str(round(covarianceParameters.A))],['parameter B ',num2str(round(covarianceParameters.B))]};

    % Plot residual gravity anomaly
    figure('Name','computeCovarianceFunction','NumberTitle','off');
    clf
    hold on
    scatter(gravityData(:,1), gravityData(:,2), 1, residuals) 
    colormap(jet)
    colorbar
    title(colorbar, 'mGal', 'FontSize', 12);
    %axis([gridParameters.MINLONG gridParameters.MAXLONG gridParameters.MINLAT gridParameters.MAXLAT])
    axis([MINLONG MAXLONG MINLAT MAXLAT])
    %ax = gca;
    %ax.PlotBoxAspectRatio = [1 abs(cos(deg2rad(mean([gridParameters.MINLAT gridParameters.MAXLAT])))) 1];
    title([gravityType,' residual gravity anomaly for block ',num2str(block_counter)])
    text((MINLONG+MAXLONG)/2,(MINLAT+MAXLAT)/2,str,'Color','black')
    xlabel('Longitude')
    ylabel('Latitude')
    %caxis([-15 15])
    saveas(gcf, [outputParameters.plotsFolder,gravityType,'Block',num2str(block_counter),'ResidualGravity.png'])

    % Plot covariance function
    figure('Name','computeCovarianceFunction','NumberTitle','off');
    clf
    hold on
    plot(rad2deg ( covarianceInfo(:,1) ), covarianceInfo(:,2), '*')
    plot(rad2deg ( covarianceInfo(:,1) ), covarianceParameters.A * fittedCovariance, '-')
    drawnow
    legend('Empirical data', 'Fitted function')
    xlabel('Spherical distance in degrees')
    ylabel('Covariance$(mGal^2)$', 'interpreter', 'latex')
    title([gravityType,' gravity auto-covariance for block ',num2str(block_counter)])    
    text(0.6,40, str,'Color','g')
    saveas(gcf,[outputParameters.plotsFolder,gravityType,'Block',num2str(block_counter),'Covariance.png'])

end