function [] = AmbiNav_Start()
%AMBINAV_START Start the AmbiNav Toolkit.
%   AMBINAV_START() first searches for and starts the 3D3A-MATLAB-Toolbox,
%   then adds all subfolders of the AmbiNav Toolkit directory to the MATLAB
%   search path.

%   ==============================================================================
%   This file is part of the 3D3A AmbiNav Toolkit.
%   
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2018 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

if exist('start3D3AMATLABToolbox','file') ~= 2
    error('Must have the 3D3A-MATLAB-Toolbox added to the search path!');
else
    start3D3AMATLABToolbox();
end

[AmbiNavDir,~,~] = fileparts(which('AmbiNav_Start'));
addpath(fullfile(AmbiNavDir, 'methods'))
addpath(fullfile(AmbiNavDir, 'rotation'))

disp('AmbiNav Toolkit found and initialized.')

end