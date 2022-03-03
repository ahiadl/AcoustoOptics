classdef Analysis
    %ANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        function this = Analysis()
        end
        
        function sigNorm = norm(sig)
            maxS = max(sig);
            minS = min(sig);
            span = maxS - minS;
            sigNorm = (sig - minS)/span;
        end
        
        function sigNorm = norm2D(sig, dims)
            if nargin == 1
                maxS = max(sig(:));
                minS = min(sig(:));
                span = maxS - minS;

                sigNorm = (sig - minS)/span;
            elseif nargin == 2             
                sizeVec = size(sig);
                elemInDims = sizeVec(dims);
                elemNonDims = sizeVec;
                elemNonDims(dims) = [];
                elemInPlanes = prod(sizeVec(dims));
                d = length(sizeVec);
                dVecOrig = 1:d;
                tmp = dVecOrig;
                tmp(dims) = [];
                dVecPerm = [dims, tmp];
                
                sigTmp = sig;
                sigTmp = permute(sigTmp, dVecPerm);
                sizeTmpVec = size(sigTmp);
                csSig = reshape(sigTmp, [elemInPlanes, sizeTmpVec(3:end) ]);
                maxS = max(csSig, [], 1);
                minS = min(csSig, [], 1);
                span = maxS - minS;
                
                repNonDims = ones(1,d-2);
                
                minMat = repmat(minS, [elemInDims, repNonDims]);
                
                sigNorm = (sigTmp - minS)/span;
                
                
            end 
        end
        
    end
end

