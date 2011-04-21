
function str = strjoin(jstr, varargin)
    str = '';
    n = size(varargin,2);
    if n >= 1, str = varargin{1}; end
    if n >= 2
        for i = 2:n
            str = [str, jstr, varargin{i}];
        end
    end
end

