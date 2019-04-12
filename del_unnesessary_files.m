

tmp=split(pwd,'/');
if ~strcmp(tmp{end},'mexInteractions')
    delete *mexa64 frmBuilder* interp1lin* intrAdsr* setclist* *slxc
else
    disp('do not delete core files')
end