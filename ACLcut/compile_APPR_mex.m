function compile_APPR_mex
% Compile APPR mex function

% different options for 32bit and 64bit Matlab
ext=mexext;
switch ext
    case {'mexw32','mexglx','mexmac','mexmaci'} %32bit
        arraydims='-compatibleArrayDims';
    case {'mexw64','mexa64','mexmaci64'} %64bit
        arraydims='-largeArrayDims';
    otherwise %potentially new architectures in the future
        warning('unknown mexext %s, assuming 64bit',ext)
        arraydims='-largeArrayDims';
end

% set up input files
location=fileparts(mfilename('fullpath'));
includes=['-I',location,'/ACLcut/APPR/matlab_matrix'];
compiles=strcat(location,'/ACLcut/APPR/',{'APPR/APPR','matlab_matrix/full','matlab_matrix/sparse'},'.cpp');
product=['/APPR.',ext];

% compile 
mex(includes,arraydims,compiles{:});

% move product to correct location
if ~strcmp(pwd,[location,'/ACLcut'])
    movefile([pwd,product],[location,'/ACLcut',product]);
end

end
