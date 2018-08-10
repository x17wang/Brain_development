function [V_def]=runFEBio_new(V,El,E,v,density,d,bcPrescribeList,bcPrescribeMagnitude,dtt)
%% FEA control settings
numTimeSteps=100; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations per time step
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retries
dtmin=dtt/3; %Minimum time step size
dtmax=dtt;  %Maximum time step size

%% CONSTRUCTING FEB MODEL
%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

% febMatID = 1:size(El,1);
% febMatID=CE;
% febMatID(CE==-2)=1;
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(fileparts(filePath)),'GIBBON-master','data','temp');


%Defining file names
febioFebFileNamePart='sphere7';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'_',num2str(d),'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'_',num2str(d),'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_',num2str(d),'_node_out.txt'];
febioLogFileName_force=[febioFebFileNamePart,'_',num2str(d),'_force_out.txt'];
febioLogFileName_vel=[febioFebFileNamePart,'_',num2str(d),'_velocity_out.txt'];

% Manually defined globals
febio_spec.Globals.Constants.T=298;
febio_spec.Globals.Constants.R=8.314e-06;
febio_spec.Globals.Constants.Fc=0;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=(1:size(El,1))'; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Tet'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(El,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=El;

%Defining node sets
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcPrescribeList;
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcPrescribeList';

%MeshData
febio_spec.MeshData.NodeData{1}.ATTR.name='nodal_loads_x';
febio_spec.MeshData.NodeData{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.MeshData.NodeData{1}.node.VAL=bcPrescribeMagnitude(:,1);
febio_spec.MeshData.NodeData{1}.node.ATTR.lid=bcPrescribeList;

febio_spec.MeshData.NodeData{2}.ATTR.name='nodal_loads_y';
febio_spec.MeshData.NodeData{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.MeshData.NodeData{2}.node.VAL=bcPrescribeMagnitude(:,2);
febio_spec.MeshData.NodeData{2}.node.ATTR.lid=bcPrescribeList;

febio_spec.MeshData.NodeData{3}.ATTR.name='nodal_loads_z';
febio_spec.MeshData.NodeData{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.MeshData.NodeData{3}.node.VAL=bcPrescribeMagnitude(:,3);
febio_spec.MeshData.NodeData{3}.node.ATTR.lid=bcPrescribeList;


%Importing initial velocity to nodes from a log file
if d>=2
%     if d < 10
%     [~, N_vel_mat,~]=importFEBio_logfile([febioFebFileName(1:end-6),'_',num2str(d-1),'_velocity_out.txt']); %Nodal displacements
    [~, N_vel_mat,~]=importFEBio_logfile(fullfile(savePath,[febioLogFileName_vel(1:end-(17+numel(num2str(d)))),num2str(d-1),'_velocity_out.txt'])); %Nodal displacements
%     end;
%     if d >= 10
%         [~, N_vel_mat,~]=importFEBio_logfile(fullfile(savePath,[febioLogFileName_vel(1:end-19),num2str(d-1),'_velocity_out.txt'])); %Nodal displacements
%     end;
    N_vel_mat=N_vel_mat(:,2:end,:);
    sizImport=size(N_vel_mat);            
    sizImport(3)=sizImport(3)+1;
    N_vel_mat_n=zeros(sizImport);
    N_vel_mat_n(:,:,2:end)=N_vel_mat;
    N_vel_mat=N_vel_mat_n;
    VN=N_vel_mat(:,:,end);
    
    febio_spec.MeshData.NodeData{4}.ATTR.name='init_vel_x';
    febio_spec.MeshData.NodeData{4}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
    febio_spec.MeshData.NodeData{4}.node.VAL=VN(:,1);
    febio_spec.MeshData.NodeData{4}.node.ATTR.lid=bcPrescribeList;
    febio_spec.MeshData.NodeData{5}.ATTR.name='init_vel_y';
    febio_spec.MeshData.NodeData{5}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
    febio_spec.MeshData.NodeData{5}.node.VAL=VN(:,2);
    febio_spec.MeshData.NodeData{5}.node.ATTR.lid=bcPrescribeList;
    febio_spec.MeshData.NodeData{6}.ATTR.name='init_vel_z';
    febio_spec.MeshData.NodeData{6}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
    febio_spec.MeshData.NodeData{6}.node.VAL=VN(:,3);
    febio_spec.MeshData.NodeData{6}.node.ATTR.lid=bcPrescribeList;

    febio_spec.Initial.init{1}.ATTR.bc='vx';
    febio_spec.Initial.init{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
    febio_spec.Initial.init{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{4}.ATTR.name;
    
    febio_spec.Initial.init{2}.ATTR.bc='vy';
    febio_spec.Initial.init{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
    febio_spec.Initial.init{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{5}.ATTR.name;
    
    febio_spec.Initial.init{3}.ATTR.bc='vz';
    febio_spec.Initial.init{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
    febio_spec.Initial.init{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{6}.ATTR.name;
end;

% % Manually defined globals
% febio_spec.Globals.Constants.ATTR.name={'T','R','Fc'};
% febio_spec.Globals.Constants.Entries={298,8.314e-06,0};

% DEFINING SPATIALLY VARYING MATERIAL SET
for i = 1:size(El,1)
    febio_spec.Material.material{i}.ATTR.type='neo-Hookean';
    febio_spec.Material.material{i}.ATTR.name=['tetra_mat_',num2str(i)];
    febio_spec.Material.material{i}.ATTR.id=i;
    febio_spec.Material.material{i}.E=E(i);
    febio_spec.Material.material{i}.v=v(i);
    febio_spec.Material.material{i}.density=density(i);
end

%Control section
febio_spec.Control.analysis.ATTR.type='dynamic';
febio_spec.Control.title='Brain analysis';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=dtt;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;

% for i = 2:5000
%     FEB_struct.Step{i+1}.Control=FEB_struct.Step{1}.Control;
% end;
%         FEB_struct.Control=FEB_struct.Step{s}.Control;

%Adding load information
febio_spec.Loads.nodal_load{1}.ATTR.bc='x';
febio_spec.Loads.nodal_load{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{1}.scale.VAL=1;
febio_spec.Loads.nodal_load{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{1}.ATTR.name;

febio_spec.Loads.nodal_load{1}.ATTR.bc='y';
febio_spec.Loads.nodal_load{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{1}.scale.VAL=1;
febio_spec.Loads.nodal_load{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{2}.ATTR.name;

febio_spec.Loads.nodal_load{1}.ATTR.bc='z';
febio_spec.Loads.nodal_load{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Loads.nodal_load{1}.ATTR.lid=bcPrescribeList;
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{1}.scale.VAL=1;
febio_spec.Loads.nodal_load{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{3}.ATTR.name;

%Load curves
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 1; dtt 1; dtt*numTimeSteps 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{3}.ATTR.file=febioLogFileName_vel;
febio_spec.Output.logfile.node_data{3}.ATTR.data='vx;vy;vz';
febio_spec.Output.logfile.node_data{3}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{3}.VAL=1:size(V,1);


%% Quick viewing of the FEBio input file structure
% febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='internal';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [~, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacements    
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    V_def=V+DN;
    
end;