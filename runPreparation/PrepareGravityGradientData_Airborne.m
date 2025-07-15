clear 

GravityGradient5DXcalibur=importdata('Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat');
Nanfilter = createNanFilter(GravityGradient5DXcalibur(:,4),470,608);%for new NSW file on 22/03/2024
GravityGradient5DXcalibur = GravityGradient5DXcalibur(~isnan(GravityGradient5DXcalibur(:,4).*Nanfilter), :);
    
GravityGradient5DOtwayMgalm=importdata('Data/GRAVITY_GRAD/OtwayMgalm.mat');
GravityGradient5DOtwayMgalm(:,5)= 1.0000e-04;
GravityGradientCombined = [GravityGradient5DXcalibur ; GravityGradient5DOtwayMgalm];
   
save('Data\processedData\XcaliburOtway.mat','GravityGradientCombined')

GravityGradient5DOtwayMgalm=importdata('Data/GRAVITY_GRAD/OtwayMgalm.mat');
GravityGradient5DOtwayMgalm(:,5)= 1.0000e-04;
save('Data\processedData\OtwayMgalm104.mat','GravityGradient5DOtwayMgalm')