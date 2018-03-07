classdef deterministicSolution
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FG;
        alpha;
        calibration;
        boundsX = struct( 'low', 0, ...
            'up', 25);
        step= 0.1;
        nSteps = 100;
        cSS;
        x;
        isNotConverged = true;
        isFixedPoint = false;
    end
    
    properties (Constant)
        maxIter = 50;
        Tolerance = 1.0e-05;
        lambda = 1; % choose lambda<1 to stabilize convergence
        zoom = 0;
    end
    
    
    
    methods
        
        function obj = deterministicSolution(calibration,cSS,boundsX,initAlpha)
            obj.calibration = calibration;
            obj.cSS = cSS;
            if nargin>2
                obj.boundsX = boundsX;
            else
                obj.boundsX.up =  (1+ obj.zoom) * obj.consNoProfit;% (1+ obj.zoom) * obj.cSS ;
                obj.boundsX.low = (1- obj.zoom) * obj.consCommitement ;
            end
            obj.step =  (obj.boundsX.up - obj.boundsX.low) /obj.nSteps;
            obj.x = gridX(obj);
            
            if nargin>3
                obj = fitFG(obj,initAlpha);
                
            else
                obj = fitFG(obj);
            end
            
            obj.isFixedPoint = checkFixedPoint(obj);
        end
        
        function mu = muCommitement(obj)
            mu  =  obj.calibration.sigma/(obj.calibration.sigma-1) ;
        end
        
        function c = consCommitement(obj)
            c  =  consNoProfit(obj)/ (muCommitement(obj)^(1/obj.calibration.gamma)) ;
        end
        
        function c = consNoProfit(obj)
            c  =  ( obj.calibration.a / (1- obj.calibration.beta* obj.calibration.theta))^(1/obj.calibration.gamma);
        end
        
        function gridx = gridX(obj,step)
            
            if nargin<2
                step =  obj.step;
            end
            
            n = ceil((obj.boundsX.up - obj.boundsX.low) / step) + 1;
            gridx  = (obj.boundsX.low + step * (0:(n-1)))';
        end
        
        function gridFull = fullGridFromAlpha(obj,alpha)
            gridFull = fullGridFromFG(obj,obj.FGfromAlpha(alpha));
        end
        
        function handleCons = consumptionPolicy(obj)
            grid = fullGrid(obj);
            handleCons = @(Cprev)   pchip(grid.x,  grid.g,  Cprev);
        end
        
        function handleMu = markupPolicy(obj)
            grid = fullGrid(obj);
            handleMu = @(Cprev)   pchip(grid.x,  grid.f,  Cprev);
        end
        
        function gridFull = fullGridFromFG(obj,FG)
            
            
            gridFull.x  = obj.x;
            gridFull.f1 = FG.f1;
            gridFull.g  = FG.g;
            %             n = size(gridFull.x,1);
            %             gridFull.f  = cumsum([FG.f0;  gridFull.f1( 1:(n-1)  ) ]) ;
            %
            %             % smoooth interpolation:
            %             gridFull.f1g = pchip(gridFull.x,  gridFull.f1,  gridFull.g);
            %             gridFull.fg  = pchip(gridFull.x,  gridFull.f ,  gridFull.g);
            [ gridFull.f,gridFull.fg,  gridFull.f1g ] = fullGridmex.fullGrid_mex( gridFull.x, FG.f0,  FG.f1 , FG.g);
            
            
            
        end
        
        function grid = fullGrid(obj)
            grid =  fullGridFromFG(obj,obj.FG);
        end
        
        function obj = fitFG(obj,alphaPrev)
            cal = obj.calibration;
            cSS = obj.cSS;
            boundsX = obj.boundsX;
            n = ceil((boundsX.up - boundsX.low) / obj.step) + 1;
            equationInputs;
            
            maxIter = obj.maxIter;
            Tolerance = obj.Tolerance;
            options = optimoptions('fminunc','Display','none');
            optionsCON = optimoptions('fmincon','Display','none','Algorithm','Interior-point','UseParallel',false); %,'algorithm','sqp'
            lambda = obj.lambda;
            isNotConverged = true;
            
            
            if nargin<2
                alphaF = zeros(n+1,1);
                alphaF(1) = obj.muCommitement;
                alphaG = cSS * ones(n,1) ;%+ cal.theta *( boundsX.low + obj.step*(1:n))'  ;
                alphaPrev = [alphaF;alphaG];
            end
            
            
            alphaF  = alphaPrev(1:(n+1));
            
            alphafLB = - 2 * ones(n+1,1);
            alphafLB(1) = 1;
            alphafUB = zeros(n+1,1);
            alphafUB(1) = obj.muCommitement;
            
            
            doubleDiagonal = [ zeros(1,n+2); -eye(n+1) zeros(n+1,1)];
            doubleDiagonal = doubleDiagonal + eye(n+2);
            gIncreasing = doubleDiagonal(2:n,1:n);
            gIncreasing = -gIncreasing;
            
            muConvex = - doubleDiagonal(3:(n+1),1:(n+1));
            
            alphaG  = alphaPrev((n+2):end);
            
            alphagLB = obj.consCommitement * ones(n,1);
            alphagUB = obj.consNoProfit * ones(n,1);
            
            h = waitbar(0,'Fitting F and G');
            
            for iter = 1:maxIter
                lossFun1 = @(alphaG) eq(2).sqLoss( obj.fullGridFromAlpha( [alphaF; alphaG ])) ;
                [alphaG1, valG,exitflag,fout] = fmincon(lossFun1,alphaG,gIncreasing,zeros(n-1,1),[],[],alphagLB,alphagUB,[],optionsCON); %fminunc(lossFun1,alphaG,options);
                
                alphaG = alphaG + lambda*(alphaG1-alphaG);
                
                lossFun2 = @(alphaF) eq(1).sqLoss( obj.fullGridFromAlpha( [alphaF; alphaG ]));
                [alphaF1, valF,exitflag,gout] =  fmincon(lossFun2,alphaF,muConvex,zeros(n-1,1),[],[],alphafLB,alphafUB,[],optionsCON); % fminunc(lossFun2,alphaF,options);
                
                alphaF = alphaF + lambda*(alphaF1-alphaF);
                
                totalLoss = (lossFun1(alphaG) + lossFun2(alphaF))/ n;
                waitbar(iter / maxIter,h,sprintf('Fitting F and G : the loss is  %1.6f',totalLoss))
                deltaARG = totalLoss; % norm(alphaPrev - [alphaF;alphaG])
                alphaPrev = [alphaF;alphaG];
                
                if deltaARG < Tolerance
                    disp(['The difference ' num2str( deltaARG)...
                        ' is below tolerance level after '...
                        num2str(iter)...
                        ' iterations' ])
                    isNotConverged = false;
                    break
                end
                
            end
            
            delete(h)
            if isNotConverged
                warning(['The difference ' num2str( deltaARG)...
                    ' cannot be reduced any further']);
                if nargin>1
                    obj = fitFG(obj);
                end
            end
            obj.alpha = [alphaF; alphaG ];
            obj.FG = obj.FGfromAlpha( [alphaF; alphaG ]) ;
            obj.isNotConverged = isNotConverged;
            
            
        end
        
        function cSS = steadyState(obj)
            grid =  fullGrid(obj);
            
            y = grid.g - grid.x;
            x = grid.x;
            
            xq = obj.gridX( obj.step/100);
            yq = pchip(x,  y,  xq);
            
            dif = abs(yq);
            [~,I]=min(dif,[],1);
            ix = I(1);
            cSS = xq(ix);
        end
        
        function isFixedPoint = checkFixedPoint(obj)
            isFixedPoint = abs(obj.steadyState - obj.cSS)<obj.step;
        end
        
        function fig = plotMarkup(obj)
            grid = obj.fullGrid;
            
            
            muSS = pchip(grid.x, grid.f,  obj.cSS);
            figure
            fig = plot(grid.x,grid.f);hold on
            line('XData', [obj.cSS obj.cSS], 'YData', [1 muSS],'Color','green');
            line('XData', [obj.consCommitement obj.cSS], 'YData', [muSS muSS],'Color','green');
            
            sigma = obj.calibration.sigma;
            muCommitement = sigma/(sigma-1);
            line('XData', [obj.consCommitement obj.boundsX.up], 'YData', [muCommitement muCommitement],'Color','red');
            line('XData', [obj.consCommitement obj.boundsX.up], 'YData', [1 1],'Color','red','LineStyle','--');
            
            
        end
        
        function fig = plotConsumption(obj)
            grid = obj.fullGrid;
            
            figure
            fig = plot(grid.x,grid.g,grid.x,grid.x,'b-.'); hold on
            line('XData', [obj.cSS obj.cSS], 'YData', [obj.boundsX.low obj.cSS ],'Color','green');
            line('XData', [obj.boundsX.low obj.cSS], 'YData', [obj.cSS  obj.cSS ],'Color','green');
            
            naiveGuess = obj.consCommitement + obj.calibration.theta*(  obj.consNoProfit - obj.consCommitement);
            
            line('XData', [obj.boundsX.low obj.boundsX.up], 'YData', [obj.consCommitement obj.consCommitement],'Color','red');
            line('XData', [obj.boundsX.low obj.boundsX.up], 'YData', [obj.consNoProfit obj.consNoProfit],'Color','red','LineStyle','--');
            
            line('XData', [naiveGuess naiveGuess], 'YData', [obj.boundsX.low naiveGuess],'Color','magenta','LineStyle','--');
            
        end
        function   disp(obj)
            iMax = size(obj,2);
            for ix = 1:iMax
                disp([~obj(ix).isNotConverged obj(ix).calibration.theta obj(ix).calibration.sigma  obj(ix).steadyState])
            end
        end
        
        
        function obj = findFixPoint(obj)
            
            
            
            
            hFP = waitbar(0,'Performing fixed point iteration ... ');
            
            
            for ix = 1 : obj.maxIter
                
                waitbar(ix / obj.maxIter,hFP)
                if obj.checkFixedPoint
                    disp('Fixed point found. ')
                    obj.isFixedPoint =  obj.checkFixedPoint
                    break
                end
                
                newGuess =  obj.cSS +  (1/ix) * (obj.steadyState - obj.cSS)
                obj = solution( obj.calibration, newGuess,obj.boundsX); %, obj.alpha
            end
            delete(hFP)
        end
        
        
        function [cGrid,cOut] =  searchSS(obj)
            grid = obj.fullGrid;
            cGrid = grid.x;
            n = size(cGrid,1);
            cal = obj.calibration;
            
            alpha = obj.alpha;
            cOut = 0*cGrid;
            
            parfor iC = 1: n
                %reverseic = n-iC+1;
                sol(iC) = solution( cal,cGrid(iC)); %,alpha
                cOut(iC) = steadyState( sol(iC) );
                
                disp([cGrid(iC) cOut(iC) sol(iC).checkFixedPoint])
                if sol(iC).checkFixedPoint
                    disp('The fixed point found')
                    %    break
                end
                % alpha = sol(iC).alpha;
                
                %  plotConsumption(sol(iC))
            end
            
            plot(cGrid,cOut,cGrid,cGrid)
        end
    end
    
    methods (Static)
        function FG =  FGfromAlpha(alpha)
            n = (size(alpha,1)-1)/2;
            FG = struct('f0',alpha(1) ,...
                'f1',alpha(2:(n+1)),...
                'g',alpha((n+2):end));
            
        end
        
        
        
        function alpha =  alphaFromFG(FG)
            n = size(FG.g,1);
            
            alpha((n+2):end) =  FG.g;
            alpha(1) = FG.f0;
            alpha(2:(n+1)) = FG.f1;
        end
        
    end
    
end

