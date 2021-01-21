% find the version of COMSOL
% currently support v41, v52, v53(a)
function COMSOL_version = findVersion(file_name)
fp = fopen(file_name,'r') ;
while (~feof(fp))
    tmp_line = fgetl(fp) ;
    if ~isempty(strfind(tmp_line,'# number of domains'))
        COMSOL_version = 41 ;
        return ;
    elseif ~isempty(strfind(tmp_line,'# number of geometric entity indices'))
        COMSOL_version = 52 ;
    elseif ~isempty(strfind(tmp_line,'# number of vertices per element'))
        COMSOL_version = 53 ;
        return;
    end
end
error('Unable to determine the version of COMSOL mphtxt file.') ;
