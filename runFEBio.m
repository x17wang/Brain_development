function [V_def]=runFEBio(V,El,E,v,density,d,modelName,bcPrescribeList,bcPrescribeMagnitude)
%% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations per time step
opt_iter=100; %Optimum number of iterations
max_retries=25; %Maximum number of retries
dtmin=0.001; %Minimum time step size
dtmax=0.01; %Maximum time step size

%% CONSTRUCTING FEB MODEL
FEB_struct.febio_spec.version='2.0';
FEB_struct.Module.Type='solid';

febMatID = 1:size(El,1);
% febMatID=CE;
% febMatID(CE==-2)=1;

% Defining file names
FEB_struct.run_filename=[modelName,'_',num2str(d),'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'_',num2str(d),'.txt']; %FEBio log file name

%Creating FEB_struct
FEB_struct.Geometry.Nodes=V;
FEB_struct.Geometry.Elements={El}; %The element sets
FEB_struct.Geometry.ElementType={'tet4'}; %The element types
FEB_struct.Geometry.ElementMat={febMatID};
FEB_struct.Geometry.ElementsPartName={'Tet'};

% Manually defined globals
FEB_struct.Globals.Constants.Names={'T','R','Fc'};
FEB_struct.Globals.Constants.Entries={298,8.314e-06,0};

% DEFINING SPATIALLY VARYING MATERIAL SET
for i = 1:size(El,1)
    FEB_struct.Materials{i}.Type='neo-Hookean';
    FEB_struct.Materials{i}.Name=['tetra_mat_',num2str(i)];
    FEB_struct.Materials{i}.Properties={'E','v'};
    FEB_struct.Materials{i}.Values={E(i),v(i)};
%             FEB_struct.Materials{i}.PropAttrName{1}='lc';
%             FEB_struct.Materials{i}.PropAttrVal{1}=i;
%             FEB_struct.Materials{i}.PropAttrVal=i;
end

%Control section
FEB_struct.Control.AnalysisType='dynamic';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Control.Values={numTimeSteps,1/numTimeSteps,...  
    max_refs,max_ups,...
    0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter'};
FEB_struct.Control.TimeStepperValues={dtmin,dtmax,max_retries,opt_iter};
%         FEB_struct.Control=FEB_struct.Step{s}.Control;

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=bcPrescribeList;
FEB_struct.Geometry.NodeSet{1}.Name='bcPrescribeList';

%Adding load information
FEB_struct.Loads.Nodal_load{1}.bc='x';
FEB_struct.Loads.Nodal_load{1}.lc=1;
% FEB_struct.Loads.Nodal_load{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Loads.Nodal_load{1}.Set=bcPrescribeList;
FEB_struct.Loads.Nodal_load{1}.nodeScale=bcPrescribeMagnitude(:,1);

FEB_struct.Loads.Nodal_load{2}.bc='y';
FEB_struct.Loads.Nodal_load{2}.lc=1;
% FEB_struct.Loads.Nodal_load{2}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Loads.Nodal_load{2}.Set=bcPrescribeList;
FEB_struct.Loads.Nodal_load{2}.nodeScale=bcPrescribeMagnitude(:,2);

FEB_struct.Loads.Nodal_load{3}.bc='z';
FEB_struct.Loads.Nodal_load{3}.lc=1;
% FEB_struct.Loads.Nodal_load{3}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Loads.Nodal_load{3}.Set=bcPrescribeList;
FEB_struct.Loads.Nodal_load{3}.nodeScale=bcPrescribeMagnitude(:,3);%Load curves

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1]};

% Load curvesdeset
% for i = 1:size(El,1)
%     FEB_struct.LoadData.LoadCurves.id(i)=i;
%     FEB_struct.LoadData.LoadCurves.type{i}='linear';
%     FEB_struct.LoadData.LoadCurves.loadPoints{i}=[0 W{i}; 1 W{i}*1.74];
% end;

%Adding output requests 
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','shell thickness'};

%Specify log file output
run_disp_output_name=[FEB_struct.run_filename(1:end-4),'_node_out.txt'];
run_force_output_name=[FEB_struct.run_filename(1:end-4),'_force_out.txt'];
FEB_struct.run_output_names={run_disp_output_name,run_force_output_name};
FEB_struct.output_types={'node_data','node_data'};
FEB_struct.data_types={'ux;uy;uz','Rx;Ry;Rz'};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars
febStruct2febFile(FEB_struct);

%% RUNNING FEBIO JOB
FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1;
FEBioRunStruct.disp_log_on=1;
FEBioRunStruct.runMode='internal';%'internal';
FEBioRunStruct.t_check=0.25; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!

if runFlag==1 %i.e. a succesful run

    %% IMPORTING NODAL DISPLACEMENT RESULTS
    %Importing nodal displacements from a log file
    [~, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements

    DN=N_disp_mat(:,2:end,end); %Final nodal displacements

    %% CREATING NODE SET IN DEFORMED STATE
    V_def=V+DN;
%     DN_magnitude=sqrt(sum(DN.^2,2));
% 
%     %Plotting the deformed model
%     [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
% 
%     hf1=cFigure;
%     title('The deformed model','FontSize',fontSize);
%     xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
% 
%     hps=patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',CF);
% 
%     view(3); axis tight;  axis equal;  grid on;
%     colormap gjet; colorbar;
%     % camlight headlight;
%     set(gca,'FontSize',fontSize);
%     drawnow;
end;
