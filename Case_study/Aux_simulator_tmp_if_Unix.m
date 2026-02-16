function [maximum_loads] = Aux_simulator_tmp_if_Unix(run_ind, wind_speed,TI,PLExp,VFlowAng,Z0,~,isUnix)
%type_of_load: k=1, Blade 1 root flapwise bending moment ; k=2, Blade 1 root edgewise bending moment
% Versatile simulator: if isUnix ==1, set for unix environment; if isUnix ==0, set for windows environment

ChangedRandomSeed = round(-2147483648 + (2147483647-(-2147483648))*rand(1)); 
    fprintf('On auxiliary simulator, for run_ind = %s, wind_speed = %7.4f,the changed RandomSeed = %s.\n', num2str(run_ind), wind_speed, num2str(ChangedRandomSeed));
 
NewRandomSeed = ChangedRandomSeed;

if isUnix ==1
    Simulator_Folder = horzcat('/tmp/', getenv('PBS_JOBID'),'/Simulator/'); %horzcat('/home/yjchoe/Workspace/Simulator', num2str(simNum),'/'); %'/home/yjchoe/Workspace/Simulator4/';     
    Copy_command ='cp ';
    Delete_command = 'rm ';
    TurbSim_executable = 'turbsim ';
    FAST_executable = 'fast_dll.exe ';
elseif isUnix ==0   %If Windows
     Simulator_Folder = horzcat(pwd,'\Imple\Simulator1\'); %'C:\Research\Implementation\Simulator4\';     
    Copy_command ='copy ';
    Delete_command = 'del ';
    TurbSim_executable = 'TurbSim.exe ';
    FAST_executable = 'FAST.exe ';   
end


%Fixed input files : These will be replicated at every run to be inputted to the same simulator
TurbSimInputFile_Fixed =horzcat(Simulator_Folder,'TurbSim_Input.inp');
FastInputFile_Fixed =horzcat(Simulator_Folder,'FAST_NRELOffshrBsline5MW.fst');
AeroDynInputFile_Fixed =horzcat(Simulator_Folder,'NRELOffshrBsline5MW_AeroDyn.ipt');

%% Temporary input files generated
TurbSimInputFile_Temp = horzcat(Simulator_Folder,'TurbSim_Input', num2str(run_ind),'.inp');
TurbSimCopyCommand = horzcat(Copy_command, TurbSimInputFile_Fixed,' ', TurbSimInputFile_Temp);
[~,~]=system(TurbSimCopyCommand);

FastInputFile_Temp = horzcat(Simulator_Folder,'FAST_NRELOffshrBsline5MW', num2str(run_ind),'.fst');
FastCopyCommand = horzcat(Copy_command, FastInputFile_Fixed,' ', FastInputFile_Temp);
[~,~]=system(FastCopyCommand);

AeroDynInputFile_Temp = horzcat(Simulator_Folder,'NRELOffshrBsline5MW_AeroDyn', num2str(run_ind),'.ipt');
AeroDynCopyCommand = horzcat(Copy_command, AeroDynInputFile_Fixed,' ', AeroDynInputFile_Temp);
[~,~]=system(AeroDynCopyCommand);

%% Temporary output files
FastOutFile_Temp = horzcat(Simulator_Folder,'FAST_NRELOffshrBsline5MW', num2str(run_ind),'.out');
TurbSimSum_Temp = horzcat(Simulator_Folder, 'TurbSim_Input', num2str(run_ind),'.sum');
TurbSimBts_Temp = horzcat(Simulator_Folder, 'TurbSim_Input', num2str(run_ind),'.bts');

%% Execution commands
TurbSimExecutionCommand =horzcat(Simulator_Folder, TurbSim_executable,TurbSimInputFile_Temp);
FastExecutionCommand = horzcat(Simulator_Folder, FAST_executable,FastInputFile_Temp);



%% Locations of lines to change
% Positions may be changed if Turbsim input file or FAST output file are changed.

% TurbSim .inp
PositionIECturbc = 32; % is the line-position to denote IECturbc in TurbSim input file.
PositionURef = 37; % is the line-position to denote URef in TurbSim input file.
PositionRandSeed1 = 4; % is the line-position to denote RandSeed1 in TurbSim input file.
PositionPLExp = 39; % is the line-position to denote PLExp in TurbSim input file.
PositionZ0 = 40; % is the line-position to denote Z0 (Surface roughness length) in TurbSim input file.
PositionVFlowAng = 26;% is the line-position to denote VFlowAng in TurbSim input file.
% AeroDyn .ipt
PositionWindFile = 10;

% FAST .fst 
PositionAERODYN =161;

% FAST .out
PositionFastData = 9; % is the line-position to skip headings of FAST output file. 


%% Change in TurbSim input file 

%Turbulence Intensity (TI) calculation
% I_ref = 0.14;
% TurbIntens = I_ref.*(0.75.*wind_speed + 5.6)./wind_speed.*100;  %I*(0.75*V +5.6)/V*100


%Write a mean wind speed at TurbSim input file 

File_TurbSimInput = fopen(TurbSimInputFile_Temp,'r+');    
for repeat = drange(1:PositionURef-1)   %to skip lines until URef
    fgetl(File_TurbSimInput);
end
fseek(File_TurbSimInput, 0, 'cof');
fprintf(File_TurbSimInput, '%7.4f', wind_speed);     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
frewind(File_TurbSimInput);


%Write a Turbulence Intensity at TurbSim input file 
for repeat = drange(1:PositionIECturbc-1)   %to skip lines until PositionIECturbc
    fgetl(File_TurbSimInput);
end
fseek(File_TurbSimInput, 0, 'cof');
fprintf(File_TurbSimInput, '%7.4f', TI*100 );     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
frewind(File_TurbSimInput);


%Write a PLExp at TurbSim input file 
for repeat = drange(1:PositionPLExp-1)   %to skip lines until PositionIECturbc
    fgetl(File_TurbSimInput);
end
fseek(File_TurbSimInput, 0, 'cof');
fprintf(File_TurbSimInput, '%7.4f', PLExp );     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
frewind(File_TurbSimInput);

% %Write a VFlowAng at Turbsim input file
for repeat = drange(1:PositionVFlowAng-1)   %to skip lines until PositionIECturbc
     fgetl(File_TurbSimInput);
end
fseek(File_TurbSimInput, 0, 'cof');
fprintf(File_TurbSimInput, '%7.4f', VFlowAng );     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
frewind(File_TurbSimInput);

%Write a Z0 (Surface roughness length [m] (or "default")) at Turbsim input file
for repeat = drange(1:PositionZ0-1)   %to skip lines until PositionIECturbc
     fgetl(File_TurbSimInput);
end
fseek(File_TurbSimInput, 0, 'cof');
fprintf(File_TurbSimInput, '%7.4f', Z0 );     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
frewind(File_TurbSimInput);

%Change the TurbSim Random seed
for repeat = drange(1:PositionRandSeed1-1)   %to skip lines until PositionRandSeed1
    fgetl(File_TurbSimInput);
end   
fseek(File_TurbSimInput, 0, 'cof');
%fprintf(File_TurbSimInput, '%s', '-2147483648');
fprintf(File_TurbSimInput, '%s', '           ');%clear the previous Random seed 
frewind(File_TurbSimInput);
for repeat = drange(1:PositionRandSeed1-1)   %to skip lines until PositionRandSeed1
    fgetl(File_TurbSimInput);
end   
fseek(File_TurbSimInput, 0, 'cof');
%NewRandomSeed = round(-2147483648 + (2147483647-(-2147483648))*rand(1));
fprintf(File_TurbSimInput, '%d', NewRandomSeed);     

% Close the connection to TurbSim input file
fclose(File_TurbSimInput);


%% Change in AeroDyn input file 
File_AeroDynInput = fopen(AeroDynInputFile_Temp,'r+');    
for repeat = drange(1:PositionWindFile-1)   %to skip lines until WindFile
    fgetl(File_AeroDynInput);
end
fseek(File_AeroDynInput, 0, 'cof');
turbsim_bts_temp_filename = horzcat('"TurbSim_Input', num2str(run_ind),'.bts"');
fprintf(File_AeroDynInput, '%s', turbsim_bts_temp_filename);     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
fclose(File_AeroDynInput);




%% Change in FAST input file 
File_FASTInput = fopen(FastInputFile_Temp,'r+');    
for repeat = drange(1:PositionAERODYN-1)   %to skip lines until WindFile
    fgetl(File_FASTInput);
end
fseek(File_FASTInput, 0, 'cof');
aerodyn_ipt_temp_filename = horzcat('"NRELOffshrBsline5MW_AeroDyn', num2str(run_ind),'.ipt"');
fprintf(File_FASTInput, '%s', aerodyn_ipt_temp_filename);     %Ref: http://www.mathworks.com/help/techdoc/ref/fprintf.html
fclose(File_FASTInput);



%% Executions
%TurbSim
[~,~]=system(TurbSimExecutionCommand);
%FAST
[~,~]=system(FastExecutionCommand);


%% Output handling

% Obtaining load variables
File_FastOutput = fopen(FastOutFile_Temp,'r');  % OR if first 60 seconds can be ignored, you may use statistics ('.sts') file generated by Crunch where maximum can be obtained
for repeat = drange(1:PositionFastData-1)   %to skip lines until data
    fgetl(File_FastOutput);
end
ncols=5; %If LE3 case, this value is 46; % number of columns = how many output variables were chosen
Fast_output = fscanf(File_FastOutput, '%f', [ncols Inf]); % NOTE: column-ordered.
%Fast_output = Fast_output'; %Transpose. However, if you want to save, deal with the un-transposed data set.

% Close the connection to FAST output file 
fclose(File_FastOutput);


%Remove first 60s of FAST output data (P. Moriarty 561T).
if size(Fast_output,2) ~= 13200 %if the output is not fully generated for any reason, run the executable on auxiliary simulator folder.
    fprintf('The FAST_output length is %d at run_ind = %s, wind_speed = %7.4f, NewRandomSeed = %s.\n', size(Fast_output,2),num2str(run_ind), wind_speed, num2str(NewRandomSeed));
    maximum_loads = Aux_simulator_tmp_if_Unix(run_ind, wind_speed,TI,PLExp,VFlowAng,Z0,NewRandomSeed,isUnix);   %Run on auxiliary simulator not the original one (running on original one repeated the same problem.)
    
else
    Fast_output = Fast_output(:,1201:13200); 


    %Calculation outside of simulation for Blade 1 root edge and flap loads:
    ip = Fast_output(3,:); %RootMxc1: Blade 1 in-plane moment (i.e., the moment caused by in-plane forces) at the blade root
    oop=Fast_output(4,:);  %RootMyc1: Blade 1 out-of-plane moment (i.e., the moment caused by out-of-plane forces) at the blade root
    BldPitch1 = Fast_output(2,:);  %BldPitch1: Blade 1 pitch angle (position)
    iptd = Fast_output(5,:); %IPDefl1: In-plane tip deflection
%    ed_c =  oop.*sin(pi()/180*BldPitch1) + ip.*cos(pi()/180*BldPitch1);    %outside calculated edgewise moment
%    fl_c =  oop.*cos(pi()/180*BldPitch1) - ip.*sin(pi()/180*BldPitch1);    %outside calculated flapwise moment
%    RootMxb1 = max(abs(ed_c));  %Blade 1 edgewise moment
%    RootMyb1 = max(abs(fl_c));  %Blade 1 flapwise moment
    
    RootMxb1 = max(abs(ip));
    RootMyb1 = max(abs(oop));
    IPDefl1 = max(abs(iptd));
    maximum_loads = [RootMyb1, RootMxb1,IPDefl1];


end
    %% Deleting temp files

     system(horzcat(Delete_command, TurbSimInputFile_Temp));    %fclose all;  %Don't use in case parallel running is interfered. %To avoid some files 'being in use' 
     system(horzcat(Delete_command, FastInputFile_Temp));
     system(horzcat(Delete_command, AeroDynInputFile_Temp));
     system(horzcat(Delete_command, FastOutFile_Temp)); 
     system(horzcat(Delete_command, TurbSimSum_Temp)); 
     system(horzcat(Delete_command, TurbSimBts_Temp));  
end