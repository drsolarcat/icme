
classdef Gsr < handle
    properties
        event;
        dht;
        axes;
        pmvab;
        branchesResidue;
        branchesLength;
        branchesCombinedResidue;
    end

    methods
        function this = Gsr(event)
            this = this@handle();
            this.init(event);
        end

        function this = init(this, event)
            this.event = event;
        end

        function this = compute(this)
            this.computeDht();
%            this.computeAxes();
%            this.computeMaps();
        end

        function this = computeDht(this)
            [this.dht.Vx, this.dht.Vy, this.dht.Vz] = gsr_loop_dht(...
                min([this.event.data.data(:).Vx])-10, 10, ...
                    max([this.event.data.data(:).Vx])+10, ...
                min([this.event.data.data(:).Vy])-10, 10, ...
                    max([this.event.data.data(:).Vy])+10, ...
                min([this.event.data.data(:).Vz])-10, 10, ...
                    max([this.event.data.data(:).Vz])+10, ...
                [this.event.data.data(:).Vx], ...
                [this.event.data.data(:).Vy], ...
                [this.event.data.data(:).Vz], ...
                [this.event.data.data(:).Bx], ...
                [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz]);

            [this.dht.Vx, this.dht.Vy, this.dht.Vz] = gsr_loop_dht(...
                this.dht.Vx-10, 1, this.dht.Vx+10, ...
                this.dht.Vy-10, 1, this.dht.Vy+10, ...
                this.dht.Vz-10, 1, this.dht.Vz+10, ...
                [this.event.data.data(:).Vx], ...
                [this.event.data.data(:).Vy], ...
                [this.event.data.data(:).Vz], ...
                [this.event.data.data(:).Bx], ...
                [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz]);

            [this.dht.Vx, this.dht.Vy, this.dht.Vz] = gsr_loop_dht(...
                this.dht.Vx-1, 0.1, this.dht.Vx+1, ...
                this.dht.Vy-1, 0.1, this.dht.Vy+1, ...
                this.dht.Vz-1, 0.1, this.dht.Vz+1, ...
                [this.event.data.data(:).Vx], ...
                [this.event.data.data(:).Vy], ...
                [this.event.data.data(:).Vz], ...
                [this.event.data.data(:).Bx], ...
                [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz]);
            this.dht.V = [this.dht.Vx, this.dht.Vy, this.dht.Vz];

            c1 = cross(-[...
                [this.event.data.data(:).Vx]; ...
                [this.event.data.data(:).Vy]; ...
                [this.event.data.data(:).Vz]], [...
                [this.event.data.data(:).Bx]; ...
                [this.event.data.data(:).By]; ...
                [this.event.data.data(:).Bz]]);

            c2 = cross(-[...
                repmat(this.dht.Vx, 1, this.event.data.count); ...
                repmat(this.dht.Vy, 1, this.event.data.count); ...
                repmat(this.dht.Vz, 1, this.event.data.count)], [...
                [this.event.data.data(:).Bx]; ...
                [this.event.data.data(:).By]; ...
                [this.event.data.data(:).Bz]]);

            cm = corrcoef([c1(1,:); c1(2,:); c1(3,:)], ...
                          [c2(1,:); c2(2,:); c2(3,:)]);

            this.dht.cc = abs(cm(1, 2));
        end

        function this = computeAxes(this)
            this.computePmva();

            [this.axes.optTheta, this.axes.optPhi] = gsr_loop_axes(...
                0, 1, 90, 0, 1, 360, ...
                this.pmvab.x, this.pmvab.y, this.pmvab.z, ...
                this.dht.V, ...
                [this.event.data.data(:).Bx], ...
                [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz], ...
                [this.event.data.data(:).Pth], ...
                this.event.config.samplingInterval);

            %{
            Vht = [this.dht.Vx; this.dht.Vy; this.dht.Vz];
            z1 = RotVecArAxe(this.pmvab.z, this.pmvab.y, ...
                this.axes.optTheta*pi/180);
            z2 = RotVecArAxe(z1, this.pmvab.z, this.axes.optPhi*pi/180);
            x2 = -(Vht-dot(Vht, z2)*z2);
            x2 = x2/norm(x2);
            y2 = cross(z2, x2);

            this.axes.x = x2;
            this.axes.y = y2;
            this.axes.z = z2;
            %}
        end

        function this = computePmva(this)
            Vht = [this.dht.Vx; this.dht.Vy; this.dht.Vz];
            Vht = Vht/norm(Vht);

            P = zeros(3, 3);
            P(1,1) = 1 - Vht(1)^2;
            P(1,2) = -Vht(1)*Vht(2);
            P(1,3) = -Vht(1)*Vht(3);
            P(2,1) = -Vht(2)*Vht(1);
            P(2,2) = 1 - Vht(2)^2;
            P(2,3) = -Vht(2)*Vht(3);
            P(3,1) = -Vht(3)*Vht(1);
            P(3,2) = -Vht(3)*Vht(2);
            P(3,3) = 1 - Vht(3)^2;

            M = getVarianceMatrix(...
                [this.event.data.data(:).Bx], ...
                [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz]);
            PMP = P*M*P;

            x1 = eigSorted(PMP, 'descend');

            this.pmvab.x = -Vht;
            this.pmvab.y = x1;
            this.pmvab.z = cross(this.pmvab.x, this.pmvab.y);
            this.pmvab.z = this.pmvab.z./norm(this.pmvab.z);
        end

        function this = loopAxis(this, minTheta, dTheta, maxTheta, ...
            minPhi, dPhi, maxPhi)

            % dHT speed vector
            Vht = [this.dht.Vx; this.dht.Vy; this.dht.Vz];

            % first guess
            thetaOpt = (minTheta+maxTheta)/2;
            phiOpt = (minPhi+maxPhi)/2;

            % angle intervals for search
            aTheta = minTheta:dTheta:maxTheta;
            aPhi = minPhi:dPhi:maxPhi;
            aThetaSize = length(aTheta);
            aPhiSize = length(aPhi);

            res = [];

            this.branchesResidue = zeros(aThetaSize, aPhiSize);

            % run search in multithreaded mode
%            move this cycle to Fortran: ./loopAxes x(1) x(2) x(3) y(1) y(2) y(3) z(1) z(2) z(3) min_theta delta_theta, max_theta, min_phi, delta_phi, max_phi
            for i = 1:aThetaSize
                theta = aTheta(i);
                z1 = RotVecArAxe(this.pmvab.z, this.pmvab.y, theta*pi/180);
                for k = 1:aPhiSize
                    phi = aPhi(k);
                    z2 = RotVecArAxe(z1, this.pmvab.z, phi*pi/180);
                    x2 = - (Vht - dot(Vht, z2)*z2);
                    x2 = x2 / norm(x2);
                    y2 = cross(z2, x2);
                    % potential and pressure in this trial coordinate system
                    curve = this.computeCurve(x2, y2, z2);
                    % calcute residue for it
                    curve.computeResidue();
                    this.branchesResidue(i,k) = curve.branchesResidue;
                    this.branchesLength(i,k) = curve.branchesLength;
                    this.combinedBranchesResidue(i,k) = curve.combinedBranchesResidue;
                end
            end

            if size(res, 1) >= 1
                % clean residue map from singularities
                if isempty(max(max(R1(R1<Inf))))
                    R1max = max(max(R2(R2<Inf)));
                else
                    R1max = max(max(R1(R1<Inf)));
                end
                R1(R1==Inf) = 1.1*R1max;
                R1(isnan(R1)) = 1.1*R1max;
                L1min = min(min(L1));
                L1max = max(max(L1));

                R2max = max(max(R2(R2<Inf)));
                R2(R2==Inf) = 1.1*R2max;
                R2(isnan(R2)) = 1.1*R2max;

                [~, index] = min(this.combinedBranchesResidue(:));
                i = floor(index/aPhiSize);
                k = index-i*aPhiSize;
                optTheta = aTheta(i);
                optPhi = aPhi(k);

                [~, i] = max(res(:, 8));
                row = res(i, 1:2);
                theta_opt = atheta(row(1)) %#ok<NOPRT>
                phi_opt = aphi(row(2)) %#ok<NOPRT>

                % another way of getting optimal angels by minimum of residue
                [~, i] = min(res(:, 3));
                %row = res(i, 1:3);
                R1min = res(i, 3);
                %theta_opt = atheta(row(1)) %#ok<NOPRT>
                %phi_opt = aphi(row(2)) %#ok<NOPRT>

                % one more way of getting optimal angels by minimum of residue in
                % wide sense
                %[~, i] = min(res(:, 5));
                %row = res(i, 1:3);
                %R1min = res(i, 3);
                %theta_opt = atheta(row(1)) %#ok<NOPRT>
                %phi_opt = aphi(row(2)) %#ok<NOPRT>

                % determine what kind of map we calculating - full or partial, the
                % partial one is calculated with finer grid
                if theta_min == 0 && theta_max == 90 && phi_min == 0 && phi_max == 360
                    ftype = 'full';
                else
                    ftype = 'partial';
                end

                [~, ymva, ~, ~] = getMinVABaxes(Bx, By, Bz);
                if dot(z, ymva) < 0
                    ymva = -ymva;
                end
                theta_mva = acos(dot(z, ymva))*180/pi;
                zxy = cross(z, ymva);
                zxy = zxy/norm(zxy);
                phi_mva = acos(dot(zxy, x))*180/pi;
                if dot(zxy, y) < 0
                    phi_mva = 360-phi_mva;
                end
                phi_mva = phi_mva-90;
                if phi_mva < 0
                    phi_mva = phi_mva+360;
                end

                [~, ymvu, ~, ~] = getMinVUBaxes(Bx, By, Bz);
                if dot(z, ymvu) < 0
                    ymvu = -ymvu;
                end
                theta_mvu = acos(dot(z, ymvu))*180/pi;
                zxy = cross(z, ymvu);
                zxy = zxy/norm(zxy);
                phi_mvu = acos(dot(zxy, x))*180/pi;
                if dot(zxy, y) < 0
                    phi_mvu = 360-phi_mvu;
                end
                phi_mvu = phi_mvu-90;
                if phi_mvu < 0
                    phi_mvu = phi_mvu+360;
                end
            end
        end

        function curve = computeCurve(this, x, y, z)
            % dHT speed vector
            Vht = [this.dht.Vx; this.dht.Vy; this.dht.Vz];

            % spatial step
            dx = -dot(Vht, x)*this.event.config.samplingInterval; % m

            % project data to selected coordinate system
            [~, By, Bz] = projection(...
                [this.event.data.data(:).Bx], [this.event.data.data(:).By], ...
                [this.event.data.data(:).Bz], x, y, z);

            A = zeros(this.event.data.count, 1);
            Pt = zeros(this.event.data.count, 1);

            % integrate the data
            for i = 1:this.event.data.count
                if i > 1 && dx > 0
                    A(i) = trapz(0:dx:dx*(i-1), By(1:i)); % Tm
                end
                Pt(i) = this.event.data.data(i).Pth+(Bz(i)^2)/2/Constant.mu; % Pa
            end

            curve = Curve(A, Pt);
        end

        function this = computeMaps(this)
        end
    end
end

