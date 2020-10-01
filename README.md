# Pharmacokinetic Modeling and Optimization of Omburtamab (8H9) Delivered Through the CSF  
## Contacts  
 Code developed by Rahul S. Yerrabelli<sup>1,\*</sup> in the lab of Nai-Kong V. Cheung, MD, PhD<sup>1</sup>.  
 <sup>1</sup>Memorial Sloan Kettering Cancer Center, New York, NY, USA.  
 <sup>\*</sup>Current address is at Carle Illinois College of Medicine, Urbana, IL, USA.  
 The modeling project spanned Jul 2017 - Oct 2020.  


## Background  
* The code here is what was used to create the analyses in our manuscript, which was accepted for publication on Sep 19, 2020 to European Journal of Nuclear Medicine and Molecular Imaging (EJNMMI). The manuscript title is "Intra-Ommaya Compartmental Radioimmunotherapy using <sup>131</sup>I-Omburtamab— Pharmacokinetic Modeling to Optimize Therapeutic Index" by Rahul S. Yerrabelli, Ping He, Edward K. Fung, Kim Kramer, Pat B. Zanzonico, John L. Humm, Hongfen Guo, Neeta Pandit-Taskar, Steven M. Larson, Nai-Kong V. Cheung.  
* Please contact us if you have questions or would like to collaborate.  


## Notes, Warnings, and Potential Sources of Confusion  
* Omburtamab was originally called 8H9. 8H9 is more commonly referred to in the code.  
* 8H9 is the engineered antibody while B7-H3 is the corresponding natural antigen on tumor cells. (Likewise for 3F8 and GD2).  
* Therapeutic index (TI) is the same therapeutic ratio (TR), which is an older term, but more frequently used in the code.  
* tumorload and R0 both represent the per-volume concentration of antigens. They only differ in their units. tumorload is mol/L, while R0 is antigens/mL. Antigen-density, NR, is related, but is in antigens/cell, which is the raw unit measured. Through calculations and estimations of cell size, you get tumorload and R0 from NR.  
* If figure files are not be created when you run the code, ensure that the appropriate folder structure is in place. The MATLAB® code will create the figure files, but will not create folders if necessary to get to the path that the figure files are supposed to be in.  


## Instructions for understanding the model MATLAB® code  
The differential equations are provided in ModOde.m. However, I recommend running solveModel.m instead, which is a wrapper function that does a multitude of extra pre-processing and post-processing steps (for example calculating AUC). This is a +Simulations folder, which continues a DefaultValues.m file, which should be run before any other code is run because it defines the default variables. Dosage.mlx file is a .mlx file I used to create the dosage response figures and are provided as a good example to learn how to run the solveModel.m and DefaultValues.m files.  

 

## Tested OS and History Details  
 Language used is MATLAB® Version 2017a on a personal Mac computer (macOS Mojave, 2018 MacBook Pro 15in). Originally I used a Windows computer when working at the MSKCC building in NYC, and then an earlier personal MacBook Pro, before this 2018 personal Mac computer. MATLAB® version might have also been different earlier.  
 Git was not used for the code files for a long time in the beginning. Additionally, the document files were not added to git until Jul 10, 2020, which is the first day I started working on preparing the manuscript for the journal EJNMMI. On that day, I also did a major restructuring of the git folder organization. After the paper was accepted for publication, the commit history and all non-essential files were stripped away before releasing the repository publicly.  
