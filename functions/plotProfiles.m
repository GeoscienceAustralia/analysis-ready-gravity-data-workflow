function plotProfiles(profileVector,dataVectorLeft,dataVectorRight,constantVariable, ...
                      profileName,dataNameLeft,dataNameRight,titleProfile,plotsFolder)
    % Input:  profileVector =
    %         dataVectorLeft =
    %         dataVectorRight =
    %         constantVariable =
    %         profileName =
    %         dataNameLeft =
    %         dataNameRight =
    %         titleProfile =
    %         plotsFolder =
    % 
    % Output: saved plot

    % Example: 
    % plotProfiles(Grav_grad(1:303,2),Grav_grad(1:303,3),Grav_grad(1:303,4),Grav_grad(1,1),'Elevation [m]','Gradient [mGal/m]','Latitude','Gradiometry',plotsFolder)
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-02.

    figure('Name','plotProfiles','NumberTitle','off');  
    colororder({'b','k'})
    yyaxis left
    plot(profileVector,dataVectorLeft,LineWidth=2)
    title([titleProfile,' along  ',num2str(constantVariable)])
    xlabel(profileName)
    ylabel(dataNameLeft)
    yyaxis right
    plot(profileVector,dataVectorRight,LineWidth=2)
    ylabel(dataNameRight)
    saveas(gcf,[plotsFolder,'plotProfiles',titleProfile,'.png'])
end