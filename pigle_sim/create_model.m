% Copyright (c) 2018, Nadav Avidor.
% All rights reserved.
% This file is part of the PIGLE - Particles Interacting in Generalized Langevin Equation simulator, subject to the 
% GNU/GPL-3.0-or-later.

% to be executed from PIGLE root folder
close_system('sl_pigle_main_current',0);
if exist('sl_pigle_main_current.slx','file')
    delete sl_pigle_main_current.slx;
end

load_system('sl_pigle_Population');
new_system('sl_pigle_main_current');
load_system('sl_pigle_main_current');
set_param('sl_pigle_main_current','PreloadFcn','load pigle_data.mat')
set_param('sl_pigle_main_current','InitFcn','load pigle_data.mat')

set_param('sl_pigle_main_current', 'StopTime', 'params.stop_time')

if params.interactions.active == 1
    add_block('simulink/Math Operations/Matrix Concatenate','sl_pigle_main_current/Matrix Concatenate')
    set_param('sl_pigle_main_current/Matrix Concatenate','NumInputs', num2str(length(params.prtcl)))
    add_block('sl_pigle_main_current/Matrix Concatenate','sl_pigle_main_current/Matrix Concatenate1')
    add_block('sl_pigle_main_current/Matrix Concatenate','sl_pigle_main_current/Matrix Concatenate2')

    add_block('simulink/Ports & Subsystems/Model','sl_pigle_main_current/Interactions')
    set_param('sl_pigle_main_current/Interactions','ModelFile','sl_interactions.slx')
    add_line('sl_pigle_main_current','Matrix Concatenate/1','Interactions/1','autorouting','smart')
    add_line('sl_pigle_main_current','Matrix Concatenate1/1','Interactions/2','autorouting','smart')
    add_line('sl_pigle_main_current','Matrix Concatenate2/1','Interactions/3','autorouting','smart')
else
    add_block('simulink/Sources/Constant','sl_pigle_main_current/Interactions');
    set_param('sl_pigle_main_current/Interactions','Value','zeros(params.model_dim,params.Nprtcl)')
    set_param('sl_pigle_main_current/Interactions','VectorParams1D','off')
end

for i=1:length(params.prtcl)
    
    iStr = num2str(i);
    
    %% 'Population block (upper-level block)
    
    add_block('sl_pigle_Population/Population','sl_pigle_main_current/Population');
    Simulink.BlockDiagram.expandSubsystem('sl_pigle_main_current/Population')

    delete_unconnected_lines(['sl_pigle_main_current'])

    if params.interactions.active == 1
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/identity ' iStr]);
        set_param(['sl_pigle_main_current/identity ' iStr],'Value',[iStr '*ones(1,params.prtcl('  iStr ').Nprtcl)'])
        set_param(['sl_pigle_main_current/identity ' iStr],'VectorParams1D','off')
        
        add_line('sl_pigle_main_current','Population 0/1',['Matrix Concatenate/' iStr],'autorouting','smart')
        add_line('sl_pigle_main_current',['identity ' iStr '/1'],['Matrix Concatenate1/' iStr],'autorouting','smart')
        add_line('sl_pigle_main_current','Population 0/4',['Matrix Concatenate2/' iStr],'autorouting','smart')
    end
    add_line('sl_pigle_main_current','Interactions/1','Selector 0/1','autorouting','smart')
    
    set_param('sl_pigle_main_current/Population 0','Name',['Population ' iStr])

    set_param('sl_pigle_main_current/ws p0','VariableName',['p' iStr])
    set_param('sl_pigle_main_current/ws p0','Name',['ws p' iStr])
    set_param('sl_pigle_main_current/ws pos0','VariableName',['pos' iStr])
    set_param('sl_pigle_main_current/ws pos0','Name',['ws pos' iStr])
    set_param('sl_pigle_main_current/ws pos_supercell_0','VariableName',['pos_supercell_' iStr])
    set_param('sl_pigle_main_current/ws pos_supercell_0','Name',['ws pos_supercell_' iStr])
    set_param(['sl_pigle_main_current/ws freeze'],'VariableName',['freeze_' iStr])
    set_param(['sl_pigle_main_current/ws freeze'],'Name',['ws freeze ' iStr])

    set_param('sl_pigle_main_current/Selector 0','Name',['Selector ' iStr])
    for j=1:i
        if j==1, tmp='1'; continue; end
        tmp = [tmp '+' 'params.prtcl(' num2str(j-1) ').Nprtcl'];
    end
    set_param(['sl_pigle_main_current/Selector ' iStr],'IndexParamArray',{'',tmp})
    set_param(['sl_pigle_main_current/Selector ' iStr],'OutputSizeArray',{'',['params.prtcl('  iStr ').Nprtcl']})

    set_param(['sl_pigle_main_current/Population ' iStr '/Data Store Memory'],'DataStoreName',['r_init_' iStr])
    set_param(['sl_pigle_main_current/Population ' iStr '/Data Store Memory'],'InitialValue',['params.prtcl('  iStr ').r_init'])
    set_param(['sl_pigle_main_current/Population ' iStr '/Data Store Read'],'DataStoreName' ,['r_init_' iStr])
    set_param(['sl_pigle_main_current/Population ' iStr '/Data Store Read1'],'DataStoreName',['r_init_' iStr])
    set_param(['sl_pigle_main_current/Population ' iStr '/Data Store Write'],'DataStoreName',['r_init_' iStr])
    
    %% sl_pigle_main_current/Population 0/r_init %
    
    set_param(['sl_pigle_main_current/Population ' iStr '/r_init/Matrix Concatenate'],'NumInputs', num2str(params.model_dim))
    set_param(['sl_pigle_main_current/Population ' iStr '/r_init/Matrix Concatenate1'],'NumInputs', num2str(params.model_dim))

    set_param(['sl_pigle_main_current/Population ' iStr '/r_init/Uniform Random Number x'],'Seed',['randi(100000,1,params.prtcl(' iStr ').Nprtcl)'])    
    set_param(['sl_pigle_main_current/Population ' iStr '/r_init/celldim x'],'Value',['repmat(params.supercell.celldim(1),1,params.prtcl(' iStr ').Nprtcl)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/r_init/Uniform Random Number y'],'Seed',['randi(100000,1,params.prtcl(' iStr ').Nprtcl)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/r_init/celldim y'],'Value',['repmat(params.supercell.celldim(2),1,params.prtcl(' iStr ').Nprtcl)'])
    if params.z_enabled > 0
        add_block(['sl_pigle_main_current/Population ' iStr '/r_init/celldim x'],['sl_pigle_main_current/Population ' iStr '/r_init/celldim z']);
        set_param(['sl_pigle_main_current/Population ' iStr '/r_init/celldim z'],'Value',['params.prtcl(' iStr ').r_init(3,:)'])
        add_line(['sl_pigle_main_current/Population ' iStr '/r_init'],'celldim z/1','Matrix Concatenate1/3','autorouting','smart')        
    end
    if params.theta_enabled > 0
        add_block('simulink/Math Operations/Product',['sl_pigle_main_current/Population ' iStr '/r_init/Product3']);        
        add_block(['sl_pigle_main_current/Population ' iStr '/r_init/celldim x'],['sl_pigle_main_current/Population ' iStr '/r_init/celldim theta']);
        add_line(['sl_pigle_main_current/Population ' iStr '/r_init'],'celldim theta/1','Product3/1','autorouting','smart')
        add_block(['sl_pigle_main_current/Population ' iStr '/r_init/Uniform Random Number x'],['sl_pigle_main_current/Population ' iStr '/r_init/Uniform Random Number theta']);
        add_line(['sl_pigle_main_current/Population ' iStr '/r_init'],'Uniform Random Number theta/1','Product3/2','autorouting','smart')
        add_line(['sl_pigle_main_current/Population ' iStr '/r_init'],'Product3/1',['Matrix Concatenate1/' num2str(params.model_dim)],'autorouting','smart')
        set_param(['sl_pigle_main_current/Population ' iStr '/r_init/celldim theta'],'Value',['repmat(params.supercell.celldim(' num2str(params.model_dim) '),1,params.prtcl(' iStr ').Nprtcl)'])
    end
    
    delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr '/r_init'])

    %% sl_pigle_main_current/Population 0/Delta R %
    if ~params.z_enabled
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/dz'])        
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/Pressure Equilibrium Reset'])
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/Selector5'])
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/Selector6'])
        delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr '/Delta R'])
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/Population ' iStr '/Delta R/reset_signal']);
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/reset_signal'],'Value',['zeros(1,params.prtcl(' iStr ').Nprtcl)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/reset_signal'],'VectorParams1D','off')
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/Population ' iStr '/Delta R/freeze_signal']);
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/freeze_signal'],'Value',['zeros(1,params.prtcl(' iStr ').Nprtcl)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/freeze_signal'],'VectorParams1D','off')
        add_line(['sl_pigle_main_current/Population ' iStr '/Delta R'],{'reset_signal/1','reset_signal/1','reset_signal/1','reset_signal/1'},{'dx/5','dy/5','dtheta/5','reset m/1'},'autorouting','smart')
        add_line(['sl_pigle_main_current/Population ' iStr '/Delta R'],{'freeze_signal/1','freeze_signal/1','freeze_signal/1','freeze_signal/1'},{'dx/6','dy/6','dtheta/6','freeze m/1'},'autorouting','smart')
        add_line(['sl_pigle_main_current/Population ' iStr '/Delta R'],'dtheta/1','Matrix Concatenate1/3','autorouting','smart')
        add_line(['sl_pigle_main_current/Population ' iStr '/Delta R'],'dtheta/2','Matrix Concatenate/3','autorouting','smart')
    else
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Pressure Equilibrium Reset/z_init'],'Value',['params.prtcl(' iStr ').r_init(3,:)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Pressure Equilibrium Reset/z_enabled'],'Value',['params.z_enabled'])
        
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dz/zeros'],'Value',['zeros(params.prtcl(' iStr ').momenta_dimension,params.prtcl(' iStr ').Nprtcl)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dz/ones'],'Value',['ones(params.prtcl(' iStr ').momenta_dimension,1)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dz/Selector'],'IndexParamArray',{['1:params.prtcl(' iStr ').momenta_dimension'],['1:params.prtcl(' iStr ').Nprtcl']})
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dz/mass'],'Value',['params.prtcl(' iStr ').mass'])
    end
    
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Matrix Concatenate'],'NumInputs', num2str(params.model_dim))
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Matrix Concatenate1'],'NumInputs', num2str(params.model_dim))
    
    if ~params.theta_enabled
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta'])        
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/Selector7'])
        delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/Band-Limited theta White Noise'])
%         delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/theta A'])
%         delete_block(['sl_pigle_main_current/Population ' iStr '/Delta R/theta B'])
        delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr '/Delta R'])
    else
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Band-Limited theta White Noise'],'Seed',['theta_seed_values_' iStr])
%         set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/theta A'],'Value',['params.prtcl(' iStr ').A_theta'])
%         set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/theta B'],'Value',['params.prtcl(' iStr ').B_theta'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Selector7'],'IndexParamArray',{num2str(params.model_dim),''})
        
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta/zeros'],'Value',['zeros(params.prtcl(' iStr ').momenta_dimension_theta,params.prtcl(' iStr ').Nprtcl)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta/ones'],'Value',['ones(params.prtcl(' iStr ').momenta_dimension_theta,1)'])
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta/Selector'],'IndexParamArray',{['1:params.prtcl(' iStr ').momenta_dimension_theta'],['1:params.prtcl(' iStr ').Nprtcl']})
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta/mass'],'Value',['params.prtcl(' iStr ').angular_mass'])
    end
    
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/Band-Limited White Noise'],'Seed',['seed_values_' iStr])    
%     set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/A'],'Value',['params.prtcl(' iStr ').A'])
%     set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/B'],'Value',['params.prtcl(' iStr ').B'])

    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dx/zeros'],'Value',['zeros(params.prtcl(' iStr ').momenta_dimension,params.prtcl(' iStr ').Nprtcl)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dx/ones'],'Value',['ones(params.prtcl(' iStr ').momenta_dimension,1)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dx/Selector'],'IndexParamArray',{['1:params.prtcl(' iStr ').momenta_dimension'],['1:params.prtcl(' iStr ').Nprtcl']})
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dx/mass'],'Value',['params.prtcl(' iStr ').mass'])

    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dy/zeros'],'Value',['zeros(params.prtcl(' iStr ').momenta_dimension,params.prtcl(' iStr ').Nprtcl)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dy/ones'],'Value',['ones(params.prtcl(' iStr ').momenta_dimension,1)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dy/Selector'],'IndexParamArray',{['1:params.prtcl(' iStr ').momenta_dimension'],['1:params.prtcl(' iStr ').Nprtcl']})
    set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dy/mass'],'Value',['params.prtcl(' iStr ').mass'])
    
    delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr '/Delta R'])

    
    %% sl_pigle_main_current/Population 0/Absolute Position %
    
    tmp='[1:2]';
    if params.theta_enabled, tmp = ['[1:2 ' num2str(params.model_dim) ']']; end
    set_param(['sl_pigle_main_current/Population ' iStr '/Absolute Position/celldim'],'Value',['repmat(params.unitcell.celldim(' tmp ')'',1,params.prtcl(' iStr ').Nprtcl)'])
    set_param(['sl_pigle_main_current/Population ' iStr '/Absolute Position/celldim1'],'Value',['repmat(params.supercell.celldim(' tmp ')'',1,params.prtcl(' iStr ').Nprtcl)'])
    
    tmp='[1:2]';
    if params.theta_enabled,  tmp = ['[1:2 ' num2str(params.model_dim) ']']; end %assuming theta is the last dim
    set_param(['sl_pigle_main_current/Population ' iStr '/Absolute Position/x,y,theta'],'IndexParamArray',{tmp,''})
    
    if ~params.z_enabled
        set_param(['sl_pigle_main_current/Population ' iStr '/Absolute Position/Matrix Concatenate'],'NumInputs','1')
        set_param(['sl_pigle_main_current/Population ' iStr '/Absolute Position/Matrix Concatenate1'],'NumInputs','1')
        delete_block(['sl_pigle_main_current/Population ' iStr '/Absolute Position/z'])
        delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr '/Absolute Position'])
    end
    
    tmp = ['[' num2str(1:params.model_dim) ']'];
    if params.z_enabled && params.theta_enabled
        tmp = '[1 2 4 3]';
    else
        
    end
    set_param(['sl_pigle_main_current/Population ' iStr '/Absolute Position/perm vec'],'Value',tmp)


    %% sl_pigle_main_current/Population 0/PES %
    target_block_name = ['sl_pigle_main_current/Population ' iStr '/PES'];
    set_param([target_block_name '/force X'],'NumberOfTableDimensions',num2str(params.model_dim))
    set_param([target_block_name '/force Y'],'NumberOfTableDimensions',num2str(params.model_dim))
    set_param([target_block_name '/force Z'],'NumberOfTableDimensions',num2str(params.model_dim))
    set_param([target_block_name '/force Theta'],'NumberOfTableDimensions',num2str(params.model_dim))

    set_param([target_block_name '/force X'],'Table',['params.prtcl(' iStr ').pes.fx'])
    set_param([target_block_name '/force Y'],'Table',['params.prtcl(' iStr ').pes.fy'])
    set_param([target_block_name '/force Z'],'Table',['params.prtcl(' iStr ').pes.fz'])
    set_param([target_block_name '/force Theta'],'Table',['params.prtcl(' iStr ').pes.ftheta'])

    set_param([target_block_name '/theta'],'IndexParamArray',{num2str(params.model_dim),''})

    if ~params.z_enabled
        delete_block([target_block_name '/force Z'])
        delete_block([target_block_name '/z'])
        delete_unconnected_lines([target_block_name ''])
        add_line([target_block_name ''],'force Theta/1','Matrix Concatenate/3')
    end
    
    set_param([target_block_name '/Matrix Concatenate'],'NumInputs', num2str(params.model_dim))
    
    if ~params.theta_enabled
        delete_block([target_block_name '/force Theta'])
        delete_block([target_block_name '/theta'])
    end
    
    delete_unconnected_lines([target_block_name ''])

    if params.theta_enabled
        % in case 'z' isn't enabled, rewrite the third breakpoint in X,Y,Theta
        set_param([target_block_name '/force Theta'],['BreakpointsForDimension' num2str(params.model_dim)],'params.unitcell.theta')
        set_param([target_block_name '/force X'],['BreakpointsForDimension' num2str(params.model_dim)],'params.unitcell.theta')
        set_param([target_block_name '/force Y'],['BreakpointsForDimension' num2str(params.model_dim)],'params.unitcell.theta')
    end
    
    if params.theta_enabled && ~params.z_enabled
        add_line([target_block_name ''],{'theta/1','theta/1','theta/1'},{'force X/3','force Y/3','force Theta/3'},'autorouting','smart')
    end

    %% sl_pigle_main_current/Population 0/Friction %
    
    target_block_name = ['sl_pigle_main_current/Population ' iStr '/Friction'];
    
    if params.prtcl(i).A_spatial_depended_friction        
        set_param([target_block_name '/A_nD_Array'],'NumberOfTableDimensions',num2str(length(size(params.prtcl(i).A))))
        set_param([target_block_name '/B_nD_Array'],'NumberOfTableDimensions',num2str(length(size(params.prtcl(i).A))))

        % An LE friction will have 1x1 dimensions, and the 'A' matrix won't
        % have extra dimensions on top of the nD matrix which is used to spatially scale
        % the friction
        if length(size(params.prtcl(i).friction.scaleMat)) == length(size(params.prtcl(i).A))
            set_param([target_block_name '/A_nD_Array'],'NumSelectionDims','0');
            set_param([target_block_name '/B_nD_Array'],'NumSelectionDims','0');
        end
        
        add_line([target_block_name ''],{'PLU_y/1','PLU_y/2','PLU_y/1','PLU_y/2'},{'A_nD_Array/1','A_nD_Array/2','B_nD_Array/1','B_nD_Array/2'},'autorouting','smart')
        add_line([target_block_name ''],{'PLU_x/1','PLU_x/2','PLU_x/1','PLU_x/2'},{'A_nD_Array/3','A_nD_Array/4','B_nD_Array/3','B_nD_Array/4'},'autorouting','smart')
        if params.z_enabled
            add_line([target_block_name ''],{'PLU_z/1','PLU_z/2','PLU_z/1','PLU_z/2'},{'A_nD_Array/5','A_nD_Array/6','B_nD_Array/5','B_nD_Array/6'},'autorouting','smart')
        end
        if params.theta_enabled
            tmpPort = 5+params.z_enabled*2; tmpPortS = num2str(tmpPort); tmpPortS1 = num2str(tmpPort+1);
            add_line([target_block_name ''],{'PLU_theta/1','PLU_theta/2','PLU_theta/1','PLU_theta/2'},{['A_nD_Array/' tmpPortS],['A_nD_Array/' tmpPortS],['B_nD_Array/' tmpPortS1],['B_nD_Array/' tmpPortS1]},'autorouting','smart')
        end

        set_param([target_block_name '/A_nD_Array'],'Table',['params.prtcl(' iStr ').A'])
        set_param([target_block_name '/B_nD_Array'],'Table',['params.prtcl(' iStr ').B'])    

        set_param([target_block_name '/theta'],'IndexParamArray',{num2str(params.model_dim),''})

        if ~params.z_enabled
            delete_block([target_block_name '/z'])
            delete_block([target_block_name '/PLU_z'])
        end

        if ~params.theta_enabled
            delete_block([target_block_name '/theta'])
            delete_block([target_block_name '/PLU_theta'])
        end

        delete_unconnected_lines([target_block_name ''])

        if params.theta_enabled && ~params.z_enabled
            add_line([target_block_name ''],{'PLU_theta/1','PLU_theta/1'},{'A_nD_Array/5','B_nD_Array/5'},'autorouting','smart')
            add_line([target_block_name ''],{'PLU_theta/2','PLU_theta/2'},{'A_nD_Array/6','B_nD_Array/6'},'autorouting','smart')
        end
        
        % Set matrix multiplication in each LE/GLE block (i.e. 'x','y','dz' in each population block) to be multiply element wise.
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dx/Matrix Multiply1'],'Multiplication','Element-wise(.*)')
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dx/Matrix Multiply2'],'Multiplication','Element-wise(.*)')
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dy/Matrix Multiply1'],'Multiplication','Element-wise(.*)')
        set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dy/Matrix Multiply2'],'Multiplication','Element-wise(.*)')
        if params.z_enabled
            set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dz/Matrix Multiply1'],'Multiplication','Element-wise(.*)')
            set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dz/Matrix Multiply2'],'Multiplication','Element-wise(.*)')
        end        
    
    else
        delete_block(target_block_name);
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/Population ' iStr '/A']);
        set_param(['sl_pigle_main_current/Population ' iStr '/A'],'Value',['params.prtcl(' iStr ').A'])
        set_param(['sl_pigle_main_current/Population ' iStr '/A'],'VectorParams1D','off')
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/Population ' iStr '/B']);
        set_param(['sl_pigle_main_current/Population ' iStr '/B'],'Value',['params.prtcl(' iStr ').B'])
        set_param(['sl_pigle_main_current/Population ' iStr '/B'],'VectorParams1D','off')
        
        delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr])
        
        add_line(['sl_pigle_main_current/Population ' iStr],'A/1','Delta R/2','autorouting','smart')
        add_line(['sl_pigle_main_current/Population ' iStr],'B/1','Delta R/3','autorouting','smart')
    end
    
    %% sl_pigle_main_current/Population 0/theta Friction %
    
    target_block_name = ['sl_pigle_main_current/Population ' iStr '/theta Friction'];

    if params.prtcl(i).A_spatial_depended_theta_friction
        set_param([target_block_name '/theta A_nD_Array'],'NumberOfTableDimensions',num2str(length(size(params.prtcl(i).A_theta))));
        set_param([target_block_name '/theta B_nD_Array'],'NumberOfTableDimensions',num2str(length(size(params.prtcl(i).A_theta))));

        % An LE friction will have 1x1 dimensions, and the 'A_theta' matrix won't
        % have extra dimensions on top of the nD matrix which is used to spatially scale
        % the ROTATIONAL friction
        if length(size(params.prtcl(i).friction.theta_scaleMat)) == length(size(params.prtcl(i).A_theta))
            set_param([target_block_name '/theta A_nD_Array'],'NumSelectionDims','0');
            set_param([target_block_name '/theta B_nD_Array'],'NumSelectionDims','0');
        end
        
        add_line([target_block_name ''],{'PLU_y/1','PLU_y/2','PLU_y/1','PLU_y/2'},{'theta A_nD_Array/1','theta A_nD_Array/2','theta B_nD_Array/1','theta B_nD_Array/2'},'autorouting','smart')
        add_line([target_block_name ''],{'PLU_x/1','PLU_x/2','PLU_x/1','PLU_x/2'},{'theta A_nD_Array/3','theta A_nD_Array/4','theta B_nD_Array/3','theta B_nD_Array/4'},'autorouting','smart')
        if params.z_enabled
            add_line([target_block_name ''],{'PLU_z/1','PLU_z/2','PLU_z/1','PLU_z/2'},{'theta A_nD_Array/5','theta A_nD_Array/6','theta B_nD_Array/5','theta B_nD_Array/6'},'autorouting','smart')
        end
        if params.theta_enabled
            tmpPort = 5+params.z_enabled*2; tmpPortS = num2str(tmpPort); tmpPortS1 = num2str(tmpPort+1);
            add_line([target_block_name ''],{'PLU_theta/1','PLU_theta/2','PLU_theta/1','PLU_theta/2'},{['theta A_nD_Array/' tmpPortS],['theta A_nD_Array/' tmpPortS],['theta B_nD_Array/' tmpPortS1],['theta B_nD_Array/' tmpPortS1]},'autorouting','smart')
        end

        set_param([target_block_name '/theta A_nD_Array'],'Table',['params.prtcl(' iStr ').A_theta'])
        set_param([target_block_name '/theta B_nD_Array'],'Table',['params.prtcl(' iStr ').B_theta'])    

        set_param([target_block_name '/theta'],'IndexParamArray',{num2str(params.model_dim),''})

        if ~params.z_enabled
            delete_block([target_block_name '/z'])
            delete_block([target_block_name '/PLU_z'])
        end

        if ~params.theta_enabled
            delete_block([target_block_name '/theta'])
            delete_block([target_block_name '/PLU_theta'])
        end

        delete_unconnected_lines([target_block_name ''])
    
        if params.theta_enabled && ~params.z_enabled
            add_line([target_block_name ''],{'PLU_theta/1','PLU_theta/1'},{'theta A_nD_Array/5','theta B_nD_Array/5'},'autorouting','smart')
            add_line([target_block_name ''],{'PLU_theta/2','PLU_theta/2'},{'theta A_nD_Array/6','theta B_nD_Array/6'},'autorouting','smart')
        end
        
        % Set matrix multiplication in each LE/GLE block (i.e. 'dtheta' in each population block) to be multiply element wise.
        if params.theta_enabled
            set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta/Matrix Multiply1'],'Multiplication','Element-wise(.*)')
            set_param(['sl_pigle_main_current/Population ' iStr '/Delta R/dtheta/Matrix Multiply2'],'Multiplication','Element-wise(.*)')
        end
        
    else
        
        delete_block(target_block_name);
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/Population ' iStr '/A_theta']);
        set_param(['sl_pigle_main_current/Population ' iStr '/A_theta'],'Value',['params.prtcl(' iStr ').A_theta'])
        set_param(['sl_pigle_main_current/Population ' iStr '/A_theta'],'VectorParams1D','off')
        add_block('simulink/Sources/Constant',['sl_pigle_main_current/Population ' iStr '/B_theta']);
        set_param(['sl_pigle_main_current/Population ' iStr '/B_theta'],'Value',['params.prtcl(' iStr ').B_theta'])
        set_param(['sl_pigle_main_current/Population ' iStr '/B_theta'],'VectorParams1D','off')
        
        delete_unconnected_lines(['sl_pigle_main_current/Population ' iStr])
        
        add_line(['sl_pigle_main_current/Population ' iStr],'A_theta/1','Delta R/4','autorouting','smart')
        add_line(['sl_pigle_main_current/Population ' iStr],'B_theta/1','Delta R/5','autorouting','smart')
    end
    
end

%if exist('AutoLayout','file'), AutoLayout('sl_pigle_main_current'); end
save_system('sl_pigle_main_current','sl_pigle_main_current')
close_system('sl_pigle_main_current',0); close_system('sl_pigle_Population',0)

addpath([pigle_path '/pigle_sim/mexInteractions/']);
convert_template_and_mex

