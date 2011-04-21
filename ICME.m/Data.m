
classdef Data < handle
    properties
        path = '';
        data;
        count;
    end

    methods
        function this = Data(pathOrData)
            this = this@handle();
            if nargin >= 1
                this.init(pathOrData);
            end
        end

        function this = init(this, pathOrData)
            if ischar(pathOrData)
                this.path = pathOrData;
            else
                this.data = pathOrData;
                this.count = length(pathOrData);
            end
        end

        function this = read(this, beginDate, endDate)
            [year, month, day, hour, minute, second, B, Bx, By, Bz, ...
            Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta] = ...
            data_read(this.path, beginDate, endDate);
            date = datenum(year, month, day, hour, minute, second);
            this.count = length(year);
            [this.data(1:this.count).date] = num2val(date);
            [this.data(:).B] = num2val(B);
            [this.data(:).Bx] = num2val(Bx);
            [this.data(:).By] = num2val(By);
            [this.data(:).Bz] = num2val(Bz);
            [this.data(:).Vp] = num2val(Vp);
            [this.data(:).Vx] = num2val(Vx);
            [this.data(:).Vy] = num2val(Vy);
            [this.data(:).Vz] = num2val(Vz);
            [this.data(:).Pth] = num2val(Pth);
            [this.data(:).Np] = num2val(Np);
            [this.data(:).Tp] = num2val(Tp);
            [this.data(:).Vth] = num2val(Vth);
            [this.data(:).beta] = num2val(beta);
        end

        function this = filter(this, beginDate, endDate)
            indices = [this.data(:).date] >= beginDate & ...
            [this.data(:).date] <= endDate;
            this.data(~indices) = [];
            this.count = length(this.data);
        end

        function that = copy(this)
            f = fieldnames(this);
            that = Data();
            for i = 1:length(f)
                that.(f{i}) = this.(f{i});
            end
        end

        function this = save(this)
        end

        function this = plot(this)
        end
    end
end

function varargout = cell2val(cellarray)
    [varargout{1:length(cellarray)}] = cellarray{:};
end

function varargout = num2val(numarray)
    cellarray = num2cell(numarray);
    [varargout{1:length(cellarray)}] = cellarray{:};
end

