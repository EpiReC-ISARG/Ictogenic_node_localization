function [CC,NN] = covm_simplified(X,Y)
% COVM generates covariance matrix
% X and Y can contain missing values encoded with NaN.
% NaN's are skipped, NaN do not result in a NaN output. 
% The output gives NaN only if there are insufficient input data
%
% COVM(X);
%      calculates the (auto-)correlation matrix of X
% COVM(X,Y);
%      calculates the crosscorrelation between X and Y
%
% Mode = 'Default' minimum or standard mode [default]
% 	C = X'*X; or X'*Y correlation matrix
%
% see also: DECOVM, XCOVF

%	$Id: covm.m 9032 2011-11-08 20:25:36Z schloegl $
%	Copyright (C) 2000-2005,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.

if nargin < 3, 
    if nargin == 1,
        Y = [];
    elseif nargin == 0
        error('Missing argument(s)');
    end
    if (~isnumeric(X) || ~isnumeric(Y))
        error('Invalid parameter(s)');
    end
else
	error('invalid input arguments');
end;

r1=size(X,1);
if ~isempty(Y)
        if r1~=size(Y,1),
                error('X and Y must have the same number of observations (rows).');
        end;
end;

if ~isempty(Y),
    CC = X'*Y;
 %   NN = real(X==X)'*real(Y==Y);
else
    CC = X'*X;
  %  tmp = real(X==X);
  %  NN  = tmp'*tmp;
end

NN = r1;
