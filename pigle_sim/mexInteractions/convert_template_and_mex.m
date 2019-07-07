
% This script copies the source files from [src_path], and compiles it at
% the current working directory, while setting the dimensions of several
% in/out ports of the 'intrAdsrbtFrcs' and 'setclist' S-functions.

 currPath = pwd;
 global pigle_path
 src_path = [pigle_path '/pigle_sim/mexInteractions/'];
 %cd(src_path);

nop=sum([params.prtcl(:).Nprtcl]);
nopStr = num2str(nop);
model_dim_Str = num2str(params.model_dim);

fTbltdDim1 = size(params.interactions.Fint,1);
fTbltdDim1Str = num2str(fTbltdDim1);

fTbltdDim2 = size(params.interactions.Fint,2);
fTbltdDim2Str = num2str(fTbltdDim2);

%% frmBuilder_sfun_intrAdsrbtFrcs

fname = 'frmBuilder_sfun_intrAdsrbtFrcs';
clear changeCellArr
%                       var  WIDTH         COL
changeCellArr(1,:)  = {'pos'     ,model_dim_Str,nopStr};
changeCellArr(2,:)  = {'x'       ,fTbltdDim1Str,'1'};
changeCellArr(3,:)  = {'fTbltd'  ,fTbltdDim1Str,fTbltdDim2Str};
changeCellArr(4,:)  = {'celldim' ,model_dim_Str,'1'};
changeCellArr(5,:)  = {'identity',nopStr       ,'1'};
changeCellArr(6,:)  = {'f_perm'  ,fTbltdDim2Str,'2'};
changeCellArr(7,:)  = {'f_func'  ,fTbltdDim2Str,'1'};
changeCellArr(8,:)  = {'freeze'  ,nopStr       ,'1'};
changeCellArr(9,:)  = {'clisti'  ,nopStr       ,'1'};
changeCellArr(10,:)  = {'clist'   ,nopStr       ,nopStr};
changeCellArr(11,:) = {'iaf'     ,model_dim_Str,nopStr};

rep_strings([fname '.c'],changeCellArr,src_path);

%% frmBuilder_sfun_intrAdsrbtFrcs_wrapper
fname = 'frmBuilder_sfun_intrAdsrbtFrcs_wrapper';
clear changeCellArr
%                       var1,    var2        addition
changeCellArr(1,:)  = {'#define','u_width',model_dim_Str};
changeCellArr(2,:)  = {'#define','y_width',nopStr};
changeCellArr(3,:)  = {'int_T','posDims[2]',[' = {' model_dim_Str ',' nopStr '};']};
changeCellArr(4,:)  = {'int_T','fTbltdDims[2]',[' = {' fTbltdDim1Str ',' fTbltdDim2Str '};']};

rep_strings1([fname '.c'],changeCellArr,src_path);

%% frmBuilder_sfun_setclist

fname = 'frmBuilder_sfun_setclist';
clear changeCellArr
%                      var  WIDTH         COL
changeCellArr(1,:) = {'pos',model_dim_Str,nopStr};
changeCellArr(2,:) = {'clisti',nopStr,'1'};
changeCellArr(3,:) = {'clist',nopStr,nopStr};

rep_strings([fname '.c'],changeCellArr,src_path);

%% frmBuilder_sfun_setclist_wrapper
fname = 'frmBuilder_sfun_setclist_wrapper';
clear changeCellArr
%                       var1,    var2        addition
changeCellArr(1,:)  = {'#define','u_width',model_dim_Str};
changeCellArr(2,:)  = {'int_T','posDims[2]',[' = {' model_dim_Str ',' nopStr '};']};

rep_strings1([fname '.c'],changeCellArr,src_path);

%%
copyfile([src_path 'frmBuilder_sfun_intrAdsrbtFrcs.tlc'],[currPath '/']);
copyfile([src_path 'intrAdsrbtFrcs.c'],[currPath '/']);
copyfile([src_path 'intrAdsrbtFrcs.h'],[currPath '/']);
copyfile([src_path 'interp1lin.c'],[currPath '/']);
copyfile([src_path 'interp1lin.h'],[currPath '/']);
mex('frmBuilder_sfun_intrAdsrbtFrcs.c','frmBuilder_sfun_intrAdsrbtFrcs_wrapper.c',[src_path 'interp1lin.c'],[src_path 'intrAdsrbtFrcs.c']);

copyfile([src_path 'frmBuilder_sfun_setclist.tlc'],[currPath '/']);
copyfile([src_path 'setclist.c'],[currPath '/']);
copyfile([src_path 'setclist.h'],[currPath '/']);
mex('frmBuilder_sfun_setclist.c','frmBuilder_sfun_setclist_wrapper.c',[src_path 'setclist.c']);

%cd(currPath);

function rep_strings(fname,changeCellArr,src_path)

fid = fopen([src_path fname]);
lines = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = lines{1};
relevant =  find(contains(lines,'#define'));
relevant1 =  find(contains(lines,'_NAME'));
relevant2 = intersect(relevant,relevant1);


for i=1:size(changeCellArr,1)
    varStr = changeCellArr{i,1};
    relevant3 =  find(contains(lines,varStr));
    ind = intersect(relevant2,relevant3);
    for j=1:length(ind)
        tmp = split(lines{ind(j)});
        if strcmp(tmp{end},varStr), ind=ind(j); break; end
    end    
    relevant2 = setdiff(relevant2,ind); % remove the found index, to not catch it again (can happen with: clist, clisti)
    ind1 = strfind(lines{ind+1},'WIDTH');
    lines{ind+1} = [lines{ind+1}(1:ind1-1) 'WIDTH' ' ' changeCellArr{i,2}];
    ind2 = strfind(lines{ind+2},'COL');
    lines{ind+2} = [lines{ind+2}(1:ind2-1) 'COL'   ' ' changeCellArr{i,3}];
end

% fname = ['tmp.' fname];
fid = fopen(fname,'w');
for i = 1:length(lines)
    fprintf(fid,'%s\n',lines{i});
end
fclose(fid);

end

function rep_strings1(fname,changeCellArr,src_path)

fid = fopen([src_path fname]);
lines = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = lines{1};

for i=1:size(changeCellArr,1)
    varStr1 = changeCellArr{i,1};
    varStr2 = changeCellArr{i,2};
    varStr3 = changeCellArr{i,3};
    relevant =  find(contains(lines,varStr1));
    relevant1 =  find(contains(lines,varStr2));
    ind = intersect(relevant,relevant1);
    lines{ind} = [varStr1 ' ' varStr2 ' ' varStr3];
end

% fname = ['tmp.' fname];
fid = fopen(fname,'w');
for i = 1:length(lines)
    fprintf(fid,'%s\n',lines{i});
end
fclose(fid);

end




