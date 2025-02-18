function datwrite(Long_outMatrix, Lat_outMatrix, resamplegeoidMatrix, OUTPUT_PARA, formattedDate)

    % Replace NaN values in the first column with 9999
    resamplegeoidMatrix(isnan(resamplegeoidMatrix)) = 999;
    
    % Make sure latitude values are positive
    Lat_outMatrix = abs(Lat_outMatrix);

    % Reshape matrices into column vectors
    Long_out_col = reshape(Long_outMatrix', [], 1);
    Lat_out_col = reshape(Lat_outMatrix', [], 1);
    resamplegeoid_col = reshape(resamplegeoidMatrix', [], 1);
    
    % Convert longitude and latitude to DM
    Long_out_col_DM = degrees2dm(Long_out_col);
    Lat_out_col_DM = degrees2dm(Lat_out_col);
   
    % Add two columns of zeros
    zero_columns = zeros(size(resamplegeoid_col, 1), 1);
    
    % Combine columns into a single matrix
    data = [resamplegeoid_col, Lat_out_col_DM, zero_columns, Long_out_col_DM, zero_columns, zero_columns, zero_columns];
    
    % Open the file for writing
    fileID = fopen([OUTPUT_PARA.Grids_name, 'AGQG_', formattedDate, '.dat'], 'w');
    
    % Write header
    fprintf(fileID, '../data/AGQG_%s                          www.ga.gov.au\n', formattedDate);
    
    % Write data in the specified format
    %for i = 1:1000
     for i = 1:size(data,1)
        fprintf(fileID, 'GEO  %.3f  S  %2d  %2.f  %.3f  E  %3d  %2.f  %.3f  %.2f  %.2f\n', ...
                data(i, 1), data(i, 2), data(i, 3), data(i, 4), data(i, 5), data(i, 6), data(i, 7), data(i, 8), data(i, 9));
    end
    
    % Close the file
    fclose(fileID);
end
