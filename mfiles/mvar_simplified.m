function [ARF,RCF,PE] = mvar_simplified(Y, Pmax)
% MVAR estimates Multi-Variate AutoRegressive model parameters
% Several estimation algorithms are implemented, all estimators 
% can handle data with missing values encoded as NaNs.  
%
% 	[AR,RC,PE] = mvar(Y, p);
%
% INPUT:
%  Y	 Multivariate data series 
%  p     Model order
%
% OUTPUT:
%  AR    multivariate autoregressive model parameter
%  RC    reflection coefficients (= -PARCOR coefficients)
%  PE    remaining error variance
%
% All input and output parameters are organized in columns, one column 
% corresponds to the parameters of one channel.
%
% Mode determines estimation algorithm. 
% Partial Correlation Estimation: Vieira-Morf [2] using unbiased covariance estimates.
% In [1] this mode was used and (incorrectly) denominated as Nutall-Strand. 
% Its the DEFAULT mode; according to [1] it provides the most accurate
% estimates.
%
% REFERENCES:
%  [1] A. Schl\"ogl, Comparison of Multivariate Autoregressive Estimators.
%       Signal processing, Elsevier B.V. (in press). 
%       available at http://dx.doi.org/doi:10.1016/j.sigpro.2005.11.007
%  [2] S.L. Marple "Digital Spectral Analysis with Applications" Prentice Hall, 1987.
%  [3] Schneider and Neumaier)
%  [4] Stijn de Waele, 2003
%
% A multivariate inverse filter can be realized with 
%   [AR,RC,PE] = mvar(Y,P);
%   e = mvfilter([eye(size(AR,1)),-AR],eye(size(AR,1)),Y);
%  
% see also: MVFILTER, MVFREQZ, COVM, SUMSKIPNAN, ARFIT2

%	$Id: mvar.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2006 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ========================================
% Edit by Janca: NaN testing was removed


% Inititialization
[N,M] = size(Y);

if nargin<2,
    Pmax=max([N,M])-1;
end;

if iscell(Y)
    Pmax = min(max(N ,M ),Pmax);
    C    = Y;
end;

[C(:,1:M),n] = covm_simplified(Y);
PE(:,1:M)  = C(:,1:M)/n;

%%%%% Partial Correlation Estimation: Vieira-Morf Method [2] with unbiased covariance estimation
%===== In [1] this algorithm is denoted "Nutall-Strand with
%unbiased covariance" =====%
F = Y;
B = Y;
PEF = PE(:,1:M);
PEB = PE(:,1:M);
for K = 1:Pmax,
    [D,n]	= covm_simplified(F(K+1:N,:),B(1:N-K,:));
    D = D/n;
    
    ARF(:,K*M+(1-M:0)) = D / PEB;
    ARB(:,K*M+(1-M:0)) = D'/ PEF;
    
    tmp        = F(K+1:N,:) - B(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
    B(1:N-K,:) = B(1:N-K,:) - F(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
    F(K+1:N,:) = tmp;
    
    for L = 1:K-1,
        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
        ARF(:,L*M+(1-M:0))   = tmp;
    end;
    
    RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
    RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
    
    [PEF,n] = covm_simplified(F(K+1:N,:),F(K+1:N,:));
    PEF = PEF/n;
    
    [PEB,n] = covm_simplified(B(1:N-K,:),B(1:N-K,:));
    PEB = PEB/n;
    
    PE(:,K*M+(1:M)) = PEF;
end