function M=spblkdiag(varargin)
% Block diagonal concatenation of matrix input arguments with sparse output.
%
%                                            |A 0 .. 0|
%   Y = SPBLKDIAG(A,B,...)  produces sparse( |0 B .. 0| )
%                                            |0 0 ..  |
%
%   Class support for inputs:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%      char, logical
%
%   See also BLKDIAG, DIAG, HORZCAT, VERTCAT

varargin=cellfun(@sparse,varargin,'UniformOutput',false);
M=blkdiag(varargin{:});


end
