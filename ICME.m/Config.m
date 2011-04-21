
classdef Config < handle
    properties
        fileFormat = '%u%s%u%u%u%u%u%u%u%u%s%16c%16c%u%d%d%d%d%u%f%f%f%u';
        dateFormat = 'yyyy-mm-dd HH:MM';
        commentStyle = 'shell';
        headerLines = 1;
        path;
        data;
        count;
    end

    properties(Constant = true)
        spacecrafts = {'STA', 'STB', 'WIND', 'ACE'};
        filenames = {'stereo_a', 'stereo_b', 'wind', 'ace'};
    end

    methods
        function this = Config(path)
            this = this@handle();
            this.init(path);
        end

        function this = init(this, path)
            this.path = path;
        end

        function this = read(this)
            [flag, quality, toPlot, toComputeGsr, toComputeMva, ...
            toComputeCylinderModel, toComputeTorusModel, ...
            toComputeHidalgoModel, toComputeOwensModel, toSave, spacecraft, ...
            beginDate, endDate, samplingInterval, minTheta, maxTheta, ...
            minPhi, maxPhi, numX, stepRatio, minY, maxY, order] = ...
            textread(this.path, this.fileFormat, ...
            'commentstyle', this.commentStyle, ...
            'headerlines', this.headerLines);
            this.count = length(flag);
            [this.data(1:this.count).flag] = num2val(flag);
            [this.data(:).quality] = quality{:};
            [this.data(:).toPlot] = num2val(toPlot);
            [this.data(:).toComputeGsr] = num2val(toComputeGsr);
            [this.data(:).toComputeMva] = (num2val(toComputeMva));
            [this.data(:).toComputeCylinderModel] = ...
            num2val(toComputeCylinderModel);
            [this.data(:).toComputeTorusModel] = num2val(toComputeTorusModel);
            [this.data(:).toComputeHidalgoModel] = ...
            num2val(toComputeHidalgoModel);
            [this.data(:).toComputeOwensModel] = num2val(toComputeOwensModel);
            [this.data(:).toSave] = num2val(toSave);
            [this.data(:).spacecraft] = spacecraft{:};
            [this.data(:).beginDate] = ...
            num2val(datenum(beginDate, this.dateFormat));
            [this.data(:).endDate] = num2val(datenum(endDate, this.dateFormat));
            [this.data(:).samplingInterval] = num2val(samplingInterval);
            [this.data(:).minTheta] = num2val(minTheta);
            [this.data(:).maxTheta] = num2val(maxTheta);
            [this.data(:).minPhi] = num2val(minPhi);
            [this.data(:).maxPhi] = num2val(maxPhi);
            [this.data(:).numX] = num2val(numX);
            [this.data(:).stepRatio] = num2val(stepRatio);
            [this.data(:).minY] = num2val(minY);
            [this.data(:).maxY] = num2val(maxY);
            [this.data(:).order] = num2val(order);

            if any([this.data.flag] == 2)
                this.data([this.data.flag] ~= 2) = [];
            else
                this.data([this.data.flag] == 0) = [];
            end
        end

        function this = filter(this, quality)
            this.data(~strcmp({this.data(:).quality}, quality)) = [];
            this.count = length(this.data);
        end
    end

    methods(Static = true)
        function filename = getDataFilename(spacecraft, samplingInterval)
            index = strmatch(spacecraft, Config.spacecrafts);
            if index
                filename = [Config.filenames{index}, '_', ...
                num2str(samplingInterval), '.dat'];
            else
                error('Unknown spacecraft');
            end
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

