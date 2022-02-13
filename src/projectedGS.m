% MIT License

% Copyright (c) 2022 Miloš Beočanin

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [x, it] = projectedGS(A, b, x0, xMin, xMax, itMax, errMax)
    rows = length(A);

    % trivial system
    if rows <= 1
        x = b/A;
        % projection
        if x < xMin
            x = xMin;
        end
        if x > xMax
            x = xMax;
        end
        return
    end

    x = x0;
    for it = 1:itMax      
        for row = 1:rows
            x(row) = 1/A(row,row)*(b(row) - A(row,1:row - 1)*x(1:row - 1) - A(row,row + 1:end)*x0(row + 1:end));
            % projection
            if x(row) < xMin(row)
                x(row) = xMin(row);
            end
            if x(row) > xMax(row)
                x(row) = xMax(row);
            end
        end
        if abs(x - x0) < errMax
            return
        end
        x0 = x;
    end
end
        
        
        
        
        
        
        
        
            
            
            
            
            
            
            
            