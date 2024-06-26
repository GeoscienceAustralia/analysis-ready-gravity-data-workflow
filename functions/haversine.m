function haversineDistance = haversine(lat1, lon1, lat2, lon2)
    % HaversineDistance computes the Haversine distance between two points on the Earth's surface.
    % Input:  lat1, lon1 = latitude and longitude of point 1 in radians
    %         lat2, lon2 = latitude and longitude of point 2 in radians
    % Output: haversineDistance = Haversine distance between the two points
    %         in radians
    % atan2 is the preferred choice than asin due to its ability to handle various cases and
    % avoid numerical instability.
    %
    % Example: see interpolateCovarianceFunction and computeSphericalEmpiricalCovariance
    %
    % Main functions
    % -
    % Written by ChatGPT3.5
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-12.

    dlat = lat1 - lat2;
    dlon = lon1 - lon2;
    % with atan2
%     a = sin(dlat / 2) .^ 2 + cos(lat1) .* cos(lat2) .* sin(dlon / 2) .^ 2;
%     haversineDistance = 2 * atan2(sqrt(a), sqrt(1 - a));
    % with asin
    sinPSIon2=sqrt(sin(dlat/2).^2+(sin(dlon/2).^2).*cos(lat1).*cos(lat2));
    haversineDistance=2*asin(sinPSIon2);

end
