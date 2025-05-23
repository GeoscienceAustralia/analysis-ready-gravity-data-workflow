function weightedFilter = createGridWeights(inputVector, DEM_PARA, GRID_PARA)
    % This function creates weights for grid blending. 
    %
    % Input:  inputVector = input vector
    %         DEM_PARA = DEM parameters: num_cols,num_rows
    %         GRID_PARA = Grid parameters: filterSize, filterRadius
    % 
    % Output: weightedFilter =  weighted filter for grid blending
    %           
    % Example: see computeLSC
    %
    % Main functions
    % - conv2
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-01. 
    
    % Make a mesh grid from the input vector
    inputGrid = reshape(inputVector, DEM_PARA.num_rows, DEM_PARA.num_cols);
    gridSize = GRID_PARA.filterSize;
    [xMesh, yMesh] = meshgrid(-gridSize:gridSize, -gridSize:gridSize);

    % Define the filter kernel
    gaussianFilter = ones(size(xMesh));
    gaussianFilter(sqrt(xMesh.^2 + yMesh.^2) > GRID_PARA.filterRadius) = 0; 
    gaussianFilter = gaussianFilter ./ sum(sum(gaussianFilter));

    % Apply the filter using convolution
    weightedFilter = conv2(inputGrid, gaussianFilter, 'same');
    weightedFilter = conv2(weightedFilter, gaussianFilter, 'same') .* inputGrid;
    weightedFilter = conv2(weightedFilter, gaussianFilter, 'same') .* inputGrid;
    weightedFilter = conv2(weightedFilter, gaussianFilter, 'same');
    weightedFilter = conv2(weightedFilter, gaussianFilter, 'same');
    weightedFilter = conv2(weightedFilter, gaussianFilter, 'same');
end
