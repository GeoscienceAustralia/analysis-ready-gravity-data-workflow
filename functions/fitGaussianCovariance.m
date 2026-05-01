function [sigma2, correlationLength, fittedCovariance] = fitGaussianCovariance( ...
    sphericalDistance, empiricalCovariance, varargin)
% fitGaussianCovariance1
% Robust Gaussian covariance fit
%
% Model:
%   C(d) = sigma2 * exp(-(d.^2)/(2*L^2))
%
% Features:
%   - bounded optimisation
%   - downweights zero-distance spike
%   - avoids range collapse
%   - stable and reproducible

    % -------------------- options
    p = inputParser;
    addParameter(p, 'anchor_sigma2', true,  @islogical);
    addParameter(p, 'min_distance', [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'do_plot', false, @islogical);
    parse(p, varargin{:});
    opt = p.Results;

    % enforce column vectors
    d = sphericalDistance(:);
    C = empiricalCovariance(:);

    % remove invalid data
    idx = isfinite(d) & isfinite(C);
    d = d(idx);
    C = C(idx);

    % sort by distance
    [d, i] = sort(d);
    C = C(i);

    % characteristic spacing
    dpos = d(d > 0);
    if isempty(dpos)
        error('Need at least one nonzero distance.');
    end
    ds = dpos(1);

    % -------------------- fit mask (avoid zero-distance spike)
    if isempty(opt.min_distance)
        dmin = ds;
    else
        dmin = opt.min_distance;
    end
    fitMask = d >= dmin;

    % -------------------- initial values
    sigma2_0 = max(C);
    if ~isfinite(sigma2_0) || sigma2_0 <= 0
        sigma2_0 = 1;
    end

    % half-height estimate for L
    half = 0.5 * sigma2_0;
    j = find(fitMask & C <= half, 1, 'first');
    if isempty(j)
        L0 = max(max(d)/8, 2*ds);
    else
        L0 = max(d(j) / sqrt(2*log(2)), 2*ds);
    end

    % bounds (important!)
    Lmin = 2*ds;
    Lmax = max(d);

    % -------------------- Gaussian model
    gaussian = @(d, s2, L) s2 .* exp(-(d.^2) ./ (2*L.^2));

    % -------------------- optimisation
    opts = optimoptions('fmincon', ...
        'Display','off', ...
        'Algorithm','interior-point', ...
        'OptimalityTolerance',1e-10, ...
        'StepTolerance',1e-12);

    if opt.anchor_sigma2
        sigma2 = sigma2_0;

        obj = @(x) sum( ...
            (gaussian(d(fitMask), sigma2, exp(x)) - C(fitMask)).^2 );

        x0 = log(L0);
        lb = log(Lmin);
        ub = log(Lmax);

        x = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);
        correlationLength = exp(x);

    else
        obj = @(x) sum( ...
            (gaussian(d(fitMask), exp(x(1)), exp(x(2))) - C(fitMask)).^2 );

        x0 = [log(sigma2_0); log(L0)];
        lb = [log(eps); log(Lmin)];
        ub = [log(10*sigma2_0); log(Lmax)];

        x = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);
        sigma2 = exp(x(1));
        correlationLength      = exp(x(2));
    end

    % fitted curve (full distance range)
    fittedCovariance = gaussian(d, sigma2, correlationLength);

    % -------------------- optional plot
    if opt.do_plot
        figure('Color','w');
        plot(rad2deg(d), C, 'k*'); hold on
        plot(rad2deg(d), fittedCovariance, 'r-', 'LineWidth',1.5)
    
        xlabel('Spherical distance (degrees)')
        ylabel('Covariance (m^2)')
        grid on
        legend('Empirical','Gaussian fit','Location','best')
        title(sprintf('\\sigma^2 = %.3g,  L(rad) = %.3g', sigma2, correlationLength))
    
        % -------- y-axis control based on anchoring
        if opt.anchor_sigma2
            % lock scale to anchored sill
            yMax = 1.05 * sigma2;
            yMin = min(C);         % preserve negative empirical values if present
        else
            % auto-scale from data
            yMax = max([C; fittedCovariance])/2;
            yMin = min([C; fittedCovariance]);
        end
    
        if yMin == yMax
            yMin = yMax - eps;
        end
        ylim([yMin yMax])
    
        ax = gca;
        ax.YAxis.Exponent = 0;
        ax.YAxis.TickLabelFormat = '%.3f';
    end

end