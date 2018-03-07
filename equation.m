classdef equation
    %equation class generates residuals
    
    properties
        formula;
        calibration;
        cSS
    end
    
    methods
        
        function obj = equation (calibration,formula,cSS)
            obj.calibration = calibration;
            obj.formula = formula;
            obj.cSS = cSS;
        end
        
        function resid = residual(obj,grid)
            resid = obj.formula(grid,obj.calibration,obj.cSS) ;
        end
        
        function sqL = sqLoss(obj,grid)
            y = residual(obj,grid);
            sqL = y'*y;
        end
    end
    
end

