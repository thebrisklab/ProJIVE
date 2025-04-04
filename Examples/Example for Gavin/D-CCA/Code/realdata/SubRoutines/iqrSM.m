function iqrange = iqrSM(data,isort) 
% IQRSM, Inter Quartile Range
%   Steve Marron's matlab function
% Inputs:
%     data  - column vector, or matrix of data, 
%                where each column is treated as a separate data set
%     isort - flag indicating need to sort:
%                0  ===>  Don't need to sort
%                        !!!  DATA ASSUMED TO BE IN INCREASING ORDER  !!!
%                1  ===>  Data unsorted, so first do a sort
%                            (default, when isort not specified)
% Output:
%     iqrange  - interquartile range(s)
%                   scalar when data is a column vector
%                   row vector (of iqr for each column) when data is a matrix
%
% Assumes path can find function:
%    cquantSM

%    Copyright (c) J. S. Marron 1996-2004


%  First decide whether or not to sort, based on number of inputs
if nargin == 1 ;
  iisort = 1 ;
          %  default is to sort, when isort unpsecified    
else ;
  iisort = isort ;
end ;

%  Do sort if needed
if iisort ~= 0 ;    %  then do a sort
  sdata = sort(data) ;
else ;
  sdata = data ;
end ;

%  now get IQR
if size(sdata,2) == 1 ;
  iqrange = cquantSM(sdata,.75,0) - cquantSM(sdata,.25,0) ;
          %  last 0 turns of need to sort
else ;
  iqrange = [] ;
  for icol = 1:size(sdata,2) ;
    iqrange = [iqrange; cquantSM(sdata(:,icol),.75,0) - ...
                          cquantSM(sdata(:,icol),.25,0)] ;
  end ;
end ;

