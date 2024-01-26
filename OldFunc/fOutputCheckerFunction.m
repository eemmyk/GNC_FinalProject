function [outputs] = fOutputCheckerFunction(xOutputfcn,optimValues,state,varargin)
    if ~isempty(xOutputfcn)
        outputs = double(abs(xOutputfcn) > 1);
        xOutputfcn

        if abs(xOutputfcn) > 1
            ang = 0
        end
    else
        outputs = 0;
    end
end

