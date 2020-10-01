# Intra-Ommaya Compartmental Radioimmunotherapy using 131I-Omburtamab— Pharmacokinetic Modeling to Optimize Therapeutic Index Code 
## Contacts
 Code developed by Rahul S. Yerrabelli* in the lab of Nai-Kong Cheung, MD, PhD.  
 Memorial Sloan Kettering Cancer Center, New York, NY, USA.  
 Jul 2017 - Oct 2020.  
 *Current location is at Carle Illinois College of Medicine, Urbana, IL, USA.


## Background
 The code here is what was used to create the analyses in our manuscript, which was accepted for publication on Sep 19, 2020 to European Journal of Nuclear Medicine and Molecular Imaging (EJNMMI).  
 Please contact us if you have questions or would like to help.  
 Note- Omburtamab was originally called 8H9.  


## Instructions for understanding the model MATLAB® code
The differential equations are provided in ModOde.m. However, I recommend running solveModel.m instead, which is a wrapper function that does a multitude of extra pre-processing and post-processing steps (for example calculating AUC). This is a +Simulations folder, which continues a DefaultValues.m file, which should be run before any other code is run because it defines the default variables. Dosage.mlx file is a .mlx file I used to create the dosage response figures and are provided as a good example to learn how to run the solveModel.m and DefaultValues.m files.  

 

## Tested OS and History Details
 Language used is MATLAB® Version 2017a on a personal Mac computer (macOS Mojave, 2018 MacBook Pro 15in). Originally I used a Windows computer when working at the MSKCC building in NYC, and then an earlier personal MacBook Pro, before this 2018 personal Mac computer. MATLAB® version might have also been different earlier.
 Git was not used for the code files for a long time in the beginning. Additionally, the document files were not added to git until Jul 10, 2020, which is the first day I started working on preparing the manuscript for the journal EJNMMI. On that day, I also did a major restructuring of the git folder organization. After the paper was accepted for publication, the commit history and all non-essential files were stripped away before releasing the repository publicly.


 
