
classdef Event < handle
    properties
        data;
        config;
        mva;
        gsr;
    end

    methods
        function this = computeGsr(this)
            gsr = Gsr(this);
            gsr.compute();
            this.gsr = gsr;
        end

        function this = computeMva(this)
            mva = Mva(this);
            mva.compute();
            this.mva = mva;
        end

        function this = computeCylinderModel(this)
        end

        function this = computeTorusModel(this)
        end

        function this = computeHidalgoModel(this)
        end

        function this = computeOwensModel(this)
        end

        function this = computeLeppingModel(this)
        end
    end
end

