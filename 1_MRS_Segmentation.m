
% set up FSL
setenv('FSLDIR', '/software/system/fsl/6.0.7.17');
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); 
setenv('PATH', [getenv('PATH') ':/software/system/fsl/6.0.7.17/share/fsl/bin']);

addpath(genpath('/software/local/spgtdjl/MATLAB/Gannet3.1-master'));
addpath /data/project/CannTeen/Scripts;

voxsegment_striatum=nan(24,4);

% start script: specifies directory MRS + T1 scans 
dirMRS =  '/data/project/CannTeen/data/Teenagers/';  
%dirMRS =  '/data/project/CannTeen/data/Adults/';  

% list of subjects to loop across 
subjects = [016002 016124 016018 016009 016004 016016 016008 016019 016003 016106 016020 016025 016007 016017 016010 016126 016011 016605 016123 016001 016022 016021 016012 016027]; 
%subjects = [026017 026124 026006 026108 026003 026125 026019 026001 026009 026122 026104 026027 026323 026026 026010 026011 026221 026018 026005 026016 026002 026007 026020 026012];

for i = 1:length(subjects) % i == index
    
    try 
        name_num = subjects(i);
        name = sprintf('%06d', name_num);  % formats the number as a 6-digit string with leading zeros
        subjdir = fullfile(dirMRS, name, 'V3');  % e.g., /.../016027/V1
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

