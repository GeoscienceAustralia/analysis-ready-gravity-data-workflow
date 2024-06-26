
# Hierarchy of functions

1. importAndFormatData.m
   - plotInputData

2. computeTerrainEffect.m
   - computeTerrainCorrection.m
      - plotTerrainCorrection.m
      - computePrismGravity.m
         - computeNagyFormula.m
   - filterDEM.m

3. computeLSC.m
   - plotOutputData.m
   - plotLevellingData.m
                   
   - computeCovarianceFunctionParameters.m
     -  plotCovarianceFunction.m
     -  computeSphericalEmpiricalCovariance.m
        - custom_grpstats.m 
        - haversine.m

     -  fitEmpiricalCovariance.m
        
   - precomputeCovarianceFunction.m   
   - interpolateCovarianceFunction.m
     - haversine.m        
   


                 
 
