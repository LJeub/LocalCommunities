ext=mexext;

switch ext
    case {'mexw32','mexglx','mexmac','mexmaci'}
        arraydims='-compatibleArrayDims';
    case {'mexw64','mexa64','mexmaci64'}
        arraydims='-largeArrayDims';
    otherwise
        warning('unknown mexext %s, assuming 64bit')
        arraydims='-largeArrayDims';
end
        

location=fileparts(mfilename('fullpath'));
includes=['-I',location,'/ACLcut/APPR/matlab_matrix'];
compile=strcat(location,'/ACLcut/APPR/',{'APPR/APPR','matlab_matrix/full','matlab_matrix/sparse'},'.cpp');
product=['/APPR.',ext];
mex(includes,arraydims,compile{:});
if ~strcmp(pwd,[location,'/ACLcut'])
    movefile([pwd,product],[location,'/ACLcut',product]);
end
