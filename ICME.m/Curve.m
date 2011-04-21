
classdef Curve < handle
    properties
        x;
        y;
        count;
        branchesResidue;
        branchLength;
        combinedBranchesResidue;
        isBranched;
        minLeftIndex;
        maxIndex;
        minRightIndex;
        f;
        df;
    end

    methods
        function this = Curve(x, y)
            this = Curve@handle();
            this.init(x, y);
        end

        function this = init(this, x, y)
            this.x = x;
            this.y = y;
            this.count = length(x);
        end

        function this = computeBranches(this)

%            % smooth transverse pressure, we smooth only pressure because
%            % magnetic potential is usually very smooth initially
%            % (integration)
%            curve = this.smooth();
%            A = curve.x;
%            Pt = curve.y;

%            % flip magnetic potential if it has global minimum
%            if abs(min(A)) > abs(max(A))
%                A = -A;
%            end

%            % index for maximum of magnetic potential
%            [~, maxIndexA] = max(A);

%            % first guess
%            isBranched = 1;
%            minLeftIndex = 1;
%            minRightIndex = this.count;
%            maxIndex = maxIndexA;

%            % if maximum is on the border then no braches
%            if maxIndexA == 1 || maxIndexA == this.count
%                isBranched = 0;
%            else
%                % first estimate for the borders
%                [minLeftA, minLeftIndexA] = min(A(1:maxIndexA));
%                [minRightA, minRightIndexA] = min(A(maxIndexA:end));
%                minRightIndexA = minRightIndexA+maxIndexA-1;

%                % align 'em
%                if minLeftA == max([minLeftA, minRightA])
%                    Ar = A(maxIndexA:end);
%                    Ar = Ar(Ar > minLeftA);
%                    [~, minIndexAr] = min(Ar);
%                    minRightIndexA = minIndexAr+maxIndexA-1;
%                else
%                    Al = A(1:maxIndexA);
%                    Al = Al(Al < minRightA);
%                    [~, maxIndexAl] = max(Al);
%                    minLeftIndexA = maxIndexAl+1;
%                end

%                % if minimum farther than maximum then no branches, hmmm...
%                % this check might not be needed
%                if minLeftIndexA >= minRightIndexA
%                    isBranched = 0;
%                else
%                    % here we'll try to find borders as they appear in
%                    % transverse pressure plot, we going in small steps
%                    % from the maximum of this curve

%                    % this is that small step
%                    delta = int32(0.05*max(size(Pt)));

%                    % searching for the maximum of the pressure curve, we
%                    % are using magnetic potential maximum as a starting
%                    % point
%                    indexPtl = max([minLeftIndexA, maxIndexA-2*delta]);
%                    indexPtr = min([minRightIndexA, maxIndexA+2*delta]);
%                    [~, maxIndexPt] = max(Pt(indexPtl:indexPtr));
%                    maxIndexPt = maxIndexA-2*delta+maxIndexPt;

%                    % if pressure maximum is somewhere outside of this test
%                    % interval then we consider such a curve not brached
%                    if (maxIndexPt == indexPtl || maxIndexPt == indexPtr) || ...
%                            (maxIndexPt >= minRightIndexA || ...
%                             maxIndexPt <= minLeftIndexA)
%                        isBranched = 0;
%                    else
%                        % get local minimums of pressure curve left to the
%                        % maximum
%                        [minLeftPt, minIndexPtl] = lmin(Pt(minLeftIndexA: ...
%                                                           maxIndexPt), 2);
%                        % iterate through local minimums
%                        minLeftIndexPt = 0;
%                        if ~isempty(minLeftPt)
%                            k = length(minLeftPt);
%                            while k >= 1
%                                index = minIndexPtl(k)+minLeftIndexA-1;
%                                index1 = max([1, index-delta]);
%                                index2 = min([this.count, index+delta]);
%                                % consider slopes on the both sides of the
%                                % local minimum
%                                leftP = polyfit(double(index1:index), ...
%                                                Pt(index1:index)', 1);
%                                rightP = polyfit(double(index:index2), ...
%                                                 Pt(index:index2)', 1);
%                                % if this local minimum is really deep then
%                                % we've found the border for the branch
%                                if leftP(1) < 0 && rightP(1) > 0
%                                    minLeftPt = minLeftPt(k);
%                                    minLeftIndexPt = index;
%                                    break
%                                end
%                                k = k-1;
%                            end
%                        end
%                        % if there were no deep minimums then use the
%                        % border detected from potential curve
%                        if minLeftIndexPt == 0
%                            minLeftPt = Pt(minLeftIndexA);
%                            minLeftIndexPt = minLeftIndexA;
%                        end

%                        % same thing to the right from the maximum
%                        [minRightPt, minIndexPtr] = lmin(Pt(maxIndexPt: ...
%                                                            minRightIndexA), 2);
%                        minRightIndexPt = 0;
%                        if ~isempty(minRightPt)
%                            k = 1;
%                            while k <= length(minRightPt);
%                                index = minIndexPtr(k)+maxIndexPt-1;
%                                leftIndex = max([1, index-delta]);
%                                rightIndex = min([this.count, index+delta]);
%                                leftP = polyfit(double(leftIndex:index), ...
%                                                       Pt(leftIndex:index)', 1);
%                                rightP = polyfit(double(index:rightIndex), ...
%                                                 Pt(index:rightIndex)', 1);
%                                if leftP(1) < 0 && rightP(1) > 0
%                                    minRightPt = minRightPt(k);
%                                    minRightIndexPt = index;
%                                    break
%                                end
%                                k = k+1;
%                            end
%                        end
%                        if minRightIndexPt == 0
%                            minRightPt = Pt(minRightIndexA);
%                            minRightIndexPt = minRightIndexA;
%                        end

%                        % align borders, so that they are on the same values of
%                        % pressure
%                        if Pt(minLeftIndexPt) >= Pt(minRightIndexPt)
%                            [~, minIndexPtr] = min(abs(Pt(maxIndexPt: ...
%                                                          minRightIndexPt)- ...
%                                                       minLeftPt));
%                            minRightIndexPt = minIndexPtr+maxIndexPt-1;
%                        else
%                            [~, minIndexPtl] = min(abs(Pt(minLeftIndexPt: ...
%                                                          maxIndexPt)- ...
%                                                       minRightPt));
%                            minLeftIndexPt = minIndexPtl+minLeftIndexPt-1;
%                        end

%                        % if we haven't even went from the maximum during our
%                        % search then there are no branches
%                        if minLeftIndexPt == maxIndexPt || ...
%                           minRightIndexPt == maxIndexPt
%                            isBranched = 0;
%                        else
%                            % align borders again a little bit
%                            if abs(A(minLeftIndexPt)) > abs(A(minRightIndexPt))
%                                d = [(1:minLeftIndexPt)', ...
%                                     floor(abs(abs(A(1:minLeftIndexPt))- ...
%                                               abs(A(minRightIndexPt))))];
%                                d = sortrows(d, [2 -1]);
%                                minLeftIndexPt = d(1, 1);
%                            else
%                                d = [(minRightIndexPt:this.count)', ...
%                                     floor(abs(abs(A(minRightIndexPt: ...
%                                                     this.count))- ...
%                                               abs(A(minLeftIndexPt))))];
%                                d = sortrows(d, [2 1]);
%                                minRightIndexPt = d(1, 1);
%                            end

%                            if minLeftIndexPt == maxIndexA || ...
%                               minRightIndexPt == maxIndexA
%                                isBranched = 0;
%                            else
%                                % finally the end! that was easy :)
%                                minLeftIndex = minLeftIndexPt;
%                                maxIndex = maxIndexA;
%                                minRightIndex = minRightIndexPt;
%                            end
%                        end
%                    end
%                end
%            end

%            this.isBranched = isBranched;
%            this.minLeftIndex = minLeftIndex;
%            this.maxIndex = maxIndex;
%            this.minRightIndex = minRightIndex;

            [this.isBranched, this.minLeftIndex, this.maxIndex, ...
                this.minRightIndex] = curve_compute_branches(this.x, this.y);
        end

        function this = computeResidue(this)
            % if there are any go on
            if this.isBranched
                % form branches
                inA = this.x(this.minLeftIndex:this.maxIndex);
                inPt = this.y(this.minLeftIndex:this.maxIndex);
                outA = this.x(this.maxIndex:this.minRightIndex);
                outPt = this.y(this.maxIndex:this.minRightIndex);

                % remove duplicates if any
                [inA, indices] = unique(inA);
                inPt = inPt(indices);
                [outA, indices] = unique(outA);
                outPt = outPt(indices);

                % length of the curve
                L = min([max(size(inA)) max(size(outA))]);

                % using big amount of points, resample to it
                Nx = 2*this.count;
                inTs = timeseries(inPt, inA, 'name', 'in');
                outTs = timeseries(outPt, outA, 'name', 'out');
                A1 = max([min(inA) min(outA)]);
                A2 = min([max(inA) max(outA)]);
                dA = (A2-A1)/(Nx-1);
                A = ((0:Nx-1)*dA)+A1;
                inTs = resample(inTs, A);
                outTs = resample(outTs, A);

                inPt = inTs.data(:);
                outPt = outTs.data(:);
                maxPt = min([max(inPt) max(outPt)]);
                minPt = max([min(inPt) min(outPt)]);

                % calculate residue
                if length(inPt) == Nx && length(outPt) == Nx
                    R = (sqrt(sum((inPt-outPt).^2)))/abs(maxPt-minPt);
                else
                    R = Inf;
                    L = 0;
                end
            else
                R = Inf;
                L = 0;
            end

            this.branchesResidue = R;
            this.branchLength = L;
            this.combinedBranchesResidue = R*Nx/2/L;

            [this.branchesResidue, this.branchLength, ...
                this.combinedBranchesResidue] =
        end

        function this = fit(this, order)
            % sort test data and construct rows of x and y
            xy = sortrows([this.x(this.minLeftIndex:this.minRightIndex), ...
                           this.y(this.minLeftIndex:this.minRightIndex)], 1);
            x = xy(:, 1);
            y = xy(:, 2);

            % amount of data points
            N = length(x);

            % polinomial fit for the data and it's derivatives
            p = polyfit(x, y, order);
            dp = polyder(p);
            ddp = polyder(dp);

            % determine the slope of the curve
            l = polyfit(x, y, 1);
            slope = sign(l(1));

            % test data for plotting, more dense than original
            delta = x(end)-x(1);
            if slope > 0
                xw = linspace(x(1)-6*delta, x(end)+delta, 100);
            else
                xw = linspace(x(1)-delta, x(end)+6*delta, 100);
            end

            % examine polinomial part
            extremums = roots(dp);
            inflections = roots(ddp);
            criticals = [extremums; inflections];
            % filter imaginary points and points outside of the data interval
            criticals = criticals(imag(criticals)==0 & criticals >= min(x) & ...
                                  criticals <= max(x));
            if ~isempty(criticals)
                xl = min(criticals);
                [~, nl] = min(abs(x-xl));
                xr = max(criticals);
                [~, nr] = min(abs(x-xr));
                switch this.config.order
                    case 2
                        if polyval(ddp, xl) > 0 % minimum
                            nm = round(mean([nl, N]));
                            if slope > 0
                                nl1 = nl;
                                nl2 = N;
                                nr1 = nm;
                                nr2 = N;
                            else
                                nl1 = 1;
                                nl2 = nm;
                                nr1 = 1;
                                nr2 = N;
                            end
                        else % maximum
                            nm = round(mean([1, nr]));
                            nl1 = 1;
                            nl2 = nm;
                            nr1 = nm;
                            nr2 = nr;
                        end
                    case 3
                        nl1 = 1;
                        nl2 = nl;
                        nr1 = nr;
                        nr2 = N;
                        % ...
                    otherwise
                        % ...
                end
            else
                nm = round(N/2);
                nl1 = 1;
                nl2 = nm;
                nr1 = nm;
                nr2 = N;
            end

            % on this stage we know all necessary intervals for curve
            % fitting

            p1l = polyfit(x(nl1:nl2), y(nl1:nl2), 1);
            dy_dx_l = p1l(1);
            p1r = polyfit(x(nr1:nr2), y(nr1:nr2), 1);
            dy_dx_r = p1r(1);
            nlm = round(mean([nl1, nl2]));
            nrm = round(mean([nr1, nr2]));
            xlm = x(nlm);
            ylm = y(nlm);
            xrm = x(nrm);
            yrm = y(nrm);
            al = dy_dx_l/ylm;
            ar = dy_dx_r/yrm;
            bl = log(ylm)-al*xlm;
            br = log(yrm)-ar*xrm;
            opts = optimset('lsqnonlin');
            opts.TolX = 1e-40;
            opts.TolFun = 1e-40;
            opts.MaxFunEvals = 500;
            opts.Display = 'none';
            if slope > 0
                fexpl = @(c, x) exp(c(1)*x+c(2));
                cexpl = [al bl];
                LBexpl = [0, -1e3];
                UBexpl = [10*abs(al), 1e3];
                fexpr = @(c, x) exp(c(1)*x+c(2))+c(3);
                cexpr = [ar br 0];
                LBexpr = [0, -1e3, -0.5];
                UBexpr = [10*abs(ar), 1e3, 0.5];
                sl = 0;
                sr = 0;
            else
                fexpl = @(c, x) exp(c(1)*x+c(2))+c(3);
                cexpl = [al bl 0];
                LBexpl = [-10*abs(al), -1e3, -0.5];
                UBexpl = [0, 1e3, 0.5];
                fexpr = @(c, x) exp(c(1)*x+c(2));
                cexpr = [ar br];
                LBexpr = [-10*abs(ar), -1e3];
                UBexpr = [0, 1e3];
                sl = 0;
                sr = 0;
            end

            % add some limits
            cexpl = lsqcurvefit(fexpl, cexpl, x(nl1:nl2+sl), y(nl1:nl2+sl), ...
                                LBexpl, UBexpl, opts);
            cexpr = lsqcurvefit(fexpr, cexpr, x(nr1+sr:nr2), y(nr1+sr:nr2), ...
                                LBexpr, UBexpr, opts);

            p1 = polyfit(x, y, 1);
            dy_dx = p1(1);
            nm = round(N/2);
            xm = x(nm);
            ym = y(nm);
            a = dy_dx/ym;
            b = log(ym)-a*xm;
            cexp = [a, b];
            fexp = @(c, x) exp(c(1)*x+c(2));
            cexp = lsqcurvefit(fexp, cexp, x, y, [], [], opts);
            dfexp = @(c, x) c(1)*exp(c(1)*x+c(2));

            xmin = min(x);
            [~, nmin] = min(abs(xw-xmin));
            xmax = max(x);
            [~, nmax] = min(abs(xw-xmax));

            w = 10/abs(xmax-xmin);
            flog = @(c, x) fexpl(cexpl, x).*logistic(x, c(1), w)+ ...
                logistic(x, c(1), -w).*polyval(p, x).*logistic(x, c(2), w)+ ...
                logistic(x, c(2), -w).*fexpr(cexpr, x);
            dflog = @(c, x) ...
                (cexpl(1)*fexpl(cexpl, x).*(1+exp(w*(x-c(1))))-w* ...
                fexpl(cexpl, x).*exp(w*(x-c(1)))).*logistic(x, c(1), w).^2+ ...
                (polyval(dp, x).*logistic(x, c(1), -w).*logistic(x, c(2), ...
                w)-polyval(p, x).*(-w*exp(-w*(x-c(1))).*(1+exp(w*(x-c(2))))+ ...
                w*exp(w*(x-c(2))).*(1+exp(-w*(x-c(1)))))).*logistic(x, c(1), ...
                -w).^2.*logistic(x, c(2), w).^2+(cexpr(1)*fexpr(cexpr, x).* ...
                (1+exp(-w*(x-c(2))))+w*fexpr(cexpr, x).*exp(-w*(x-c(2)))).* ...
                logistic(x, c(2), -w).^2;
            LB = [xmin-2, xmax-2];
            UB = [xmin+2, xmax+2];
            clog = [xmin, xmax];
            clog = lsqcurvefit(flog, clog, xw, ...
                [fexpl(cexpl, xw(1:nmin)), polyval(p, xw(nmin+1:nmax)), ...
                 fexpr(cexpr, xw(nmax+1:end))], LB, UB, opts);

            flog = @(x) flog(clog, x);
            dflog = @(x) dflog(clog, x);
            fexp = @(x) fexp(cexp, x);
            dfexp = @(x) dfexp(cexp, x);

            if mean(abs(flog(this.x)-this.y).^2) < ...
                mean(abs(fexp(this.x)-this.y).^2)
                this.f = flog;
                this.df = dflog;
            else
                this.f = fexp;
                this.df = dfexp;
            end
        end

        function curve = smooth(this)
            % smoothing, makes search for invariant axis a little bit
            % easier
            x = smooth(this.x, 3, 'sgolay');
            y = smooth(this.y, 3, 'sgolay');
            y = smooth(y, 0.07*this.count);
            y = smooth(y, 0.1*this.count); % Pa
            curve = Curve(x, y);
        end

        function this = plot(this)
        end
    end
end

function H = logistic(x, x0, t)
    H = 1./(1+exp(t*(x-x0)));
end

