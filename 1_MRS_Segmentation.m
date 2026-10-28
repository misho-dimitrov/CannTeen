
% set up FSL
setenv('FSLDIR', '/software/system/fsl/6.0.7.17');
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); 
setenv('PATH', [getenv('PATH') ':/software/system/fsl/6.0.7.17/share/fsl/bin']);

addpath(genpath('/software/local/<username>/MATLAB/Gannet3.1-master'));
addpath /data/project/CannTeen/Scripts;

voxsegment_striatum=nan(24,4);

% start script: specifies directory MRS + T1 scans 
dirMRS =  '/data/project/CannTeen/data/Teenagers/';  
%dirMRS =  '/data/project/CannTeen/data/Adults/';  

% list of subjects to loop across
% NOTE: actual participant IDs redacted for the public repo; placeholders shown below, replace with the real subject ID list locally
subjects = [000001 000002 000003 000004 000005 000006 000007 000008 000009 000010 000011 000012 000013 000014 000015 000016 000017 000018 000019 000020 000021 000022 000023 000024]; 
%subjects = [000001 000002 000003 000004 000005 000006 000007 000008 000009 000010 000011 000012 000013 000014 000015 000016 000017 000018 000019 000020 000021 000022 000023 000024];

for i = 1:length(subjects) % i == index
    
    try 
        name_num = subjects(i);
        name = sprintf('%06d', name_num);  % formats the number as a 6-digit string with leading zeros
        subjdir = fullfile(dirMRS, name, 'V3');  % e.g., /.../000001/V1
        cd(subjdir); % move into subject dir     
     
                %% Main function you want to run:     
             %% Striatum      
                [gm wm csf] = GannetSegmentSiemens('metab.rda', 'structural.nii'); 
                % colect all the data across participants
                voxsegment_striatum(i,1) = name_num;
                voxsegment_striatum(i,2) = gm;
                voxsegment_striatum(i,3) = wm;
                voxsegment_striatum(i,4) = csf;    
   catch err
        fprintf('Error processing subject %s: %s\n', num2str(name), err.message);
   end
    
end

cd(dirMRS);


save CannTeen_MRS_SEGMENT_V3 voxsegment_striatum

