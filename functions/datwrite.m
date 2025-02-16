function datwrite(Long_outMatrix, Lat_outMatrix, resamplegeoidMatrix, OUTPUT_PARA, formattedDate)

    % Function to convert decimal degrees to DMS
    convertToDMS = @(deg) [fix(deg), fix((deg - fix(deg)) * 60), ((deg - fix(deg)) * 60 - fix((deg - fix(deg)) * 60)) * 60];
    
    % Replace NaN values in the first column with 9999
    resamplegeoidMatrix(isnan(resamplegeoidMatrix)) = 999;
    
    % Make sure latitude values are positive
    Lat_outMatrix = abs(Lat_outMatrix);
    
    % Convert longitude and latitude to DMS
    Long_out_DMS = arrayfun(convertToDMS, Long_outMatrix, 'UniformOutput', false);
    Lat_out_DMS = arrayfun(convertToDMS, Lat_outMatrix, 'UniformOutput', false);
    
    % Reshape matrices into column vectors
    Long_out_DMS_col = cell2mat(reshape(Long_out_DMS, [], 1));
    Lat_out_DMS_col = cell2mat(reshape(Lat_out_DMS, [], 1));
    resamplegeoid_col = reshape(resamplegeoidMatrix, [], 1);
    
    % Add two columns of zeros
    zero_columns = zeros(size(resamplegeoid_col, 1), 2);
    
    % Combine columns into a single matrix
    data = [resamplegeoid_col, Lat_out_DMS_col, Long_out_DMS_col, zero_columns];
    
    % Open the file for writing
    fileID = fopen([OUTPUT_PARA.Grids_name, 'AGQG_', formattedDate, '.dat'], 'w');
    
    % Write header
    fprintf(fileID, '../data/AGQG_%s                          www.ga.gov.au\n', formattedDate);
    
    % Write data in the specified format
    for i = 1:1000
        fprintf(fileID, 'GEO  %.3f  S  %2d  %2d  %.3f  E  %3d  %2d  %.3f  %.2f  %.2f\n', ...
                data(i, 1), data(i, 2), data(i, 3), data(i, 4), data(i, 5), data(i, 6), data(i, 7), data(i, 8), data(i, 9));
    end
    
    % Close the file
    fclose(fileID);
end
