classdef Exoskeleton < handle
    properties(Access = public)
        k;
    end
    methods(Access = public)
        function e = Exoskeleton(k)
            %normalize spring constant
            e.k = k/FootDropModel.tibialisLength(pi/2);
        end
        
        function result = force(e, lM)
            result = zeros(size(lM));
            for i = 1:length(lM)
                if(lM(i) < 1)
                    result(i) = 0;    
                else
                    result(i) = e.k*(lM(i)-1);
                end
            end
        end
    end
end