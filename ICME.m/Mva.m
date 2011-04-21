
classdef Mva < handle
    properties
        event;
        mvab;
        mvub;
    end

    methods
        function this = Mva(event)
            this = this@handle();
            this.init(event);
        end

        function this = init(this, event)
            this.event = event;
        end

        function this = compute(this)
            this.mvab = Mva.computeStruct(...
                [this.event.data.data(:).Bx], ...
                [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz]);
            this.mvub = Mva.computeStruct(...
                [this.event.data.data(:).Bx]./[this.event.data.data(:).B], ...
                [this.event.data.data(:).By]./[this.event.data.data(:).B], ...
                [this.event.data.data(:).Bz]./[this.event.data.data(:).B]);
        end
    end

    methods(Static = true)
        function mvaStruct = computeStruct(Bx, By, Bz)
            mvaStruct.matrix = getVarianceMatrix(Bx, By, Bz);
            [mvaStruct.x, mvaStruct.y, mvaStruct.z, ~, lambda2, lambda3] = ...
                eigSorted(mvaStruct.matrix, 'descend');
            mvaStruct.criterion = lambda2/lambda3;
            T = mean(Bx)*mvaStruct.y(1)+ ...
                mean(By)*mvaStruct.y(1)+ ...
                mean(Bz)*mvaStruct.y(1);
            if T < 0
                mvaStruct.x = -mvaStruct.x;
                mvaStruct.y = -mvaStruct.y;
                mvaStruct.z = -mvaStruct.z;
            end
        end
    end
end

