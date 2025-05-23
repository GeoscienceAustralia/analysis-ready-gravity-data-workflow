function Nanfilter = createNanFilter(GradientVector,num_cols,num_rows)
    % Creates a nan filter based on input data
    %
    % Input:    GradientVector = Input data (vector)
    % Output:   Nanfilter = Nan filter vector
    %

    % Initialize Nan filter
    Nanfilter = 1 + zeros(size(GradientVector));

    % Set NaN values to NaN in filter
    Nanfilter(isnan(GradientVector)) = NaN;

    % Reshape Nan filter
    Nanmask_mat = reshape(Nanfilter, num_cols,num_rows);
    
    % Set NaN values in first column to NaN
    Nanmask_mat(:, 1) = NaN;

    % Define convolution kernel
    gf = ones(6, 6);

    % Apply convolution to Nan filter
    Nanmask_mat = conv2(Nanmask_mat, gf, 'same');

    % Set non-NaN values in convolution result to 1
    Nanmask_mat(~isnan(Nanmask_mat)) = 1;

    % Reshape Nan filter back to vector
    Nanfilter = Nanmask_mat(:);
end
