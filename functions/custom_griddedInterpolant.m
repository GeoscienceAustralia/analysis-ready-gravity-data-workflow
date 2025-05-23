classdef custom_griddedInterpolant
    properties
        Grids
        Values
        Method
    end

    methods
        function obj = custom_griddedInterpolant(grids, values, method)
            % Constructor to initialize the grids and values
            if nargin < 3
                method = 'linear';
            end
            if ~iscell(grids)
                error('Grids should be a cell array.');
            end
            if ~isnumeric(values)
                error('Values should be numeric.');
            end
            obj.Grids = grids;
            obj.Values = values;
            obj.Method = method;
        end

        function Vq = interp(obj, varargin)
            % Perform interpolation using interpn
            if length(varargin) ~= length(obj.Grids)
                error('Number of query grids must match the number of original grids.');
            end
            Vq = interpn(obj.Grids{:}, obj.Values, varargin{:}, obj.Method);
        end
    end
end
