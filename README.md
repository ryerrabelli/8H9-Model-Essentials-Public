# Pharmacokinetic Modeling and Optimization of Omburtamab (8H9) Delivered Through the CSF  

**NOTE: For the most updated version of this document and all associated code and files, see: https://github.com/ryerrabelli/8H9-Model-Essentials-Public**

## Contact  

 Developed by [Rahul S. Yerrabelli](https://orcid.org/0000-0002-7670-9601)<sup>1,2</sup> in the lab of [Nai-Kong V. Cheung, MD, PhD](https://orcid.org/0000-0001-6323-5171
)<sup>1</sup>.  
 1. [Memorial Sloan Kettering Cancer Center, New York, NY, USA](https://www.mskcc.org/research-areas/labs/nai-kong-cheung)  
 1. [Carle Illinois College of Medicine, University of Illinois at Urbana-Champaign, Urbana, IL, USA](https://medicine.illinois.edu/).  
 


## Background  
* The code here is what I used to create the analyses in our manuscript, which was accepted for publication on Sep 19, 2020 to European Journal of Nuclear Medicine and Molecular Imaging (EJNMMI). The manuscript title is ["Intra-Ommaya Compartmental Radioimmunotherapy using <sup>131</sup>I-Omburtamab— Pharmacokinetic Modeling to Optimize Therapeutic Index" by Rahul S. Yerrabelli, Ping He, Edward K. Fung, Kim Kramer, Pat B. Zanzonico, John L. Humm, Hongfen Guo, Neeta Pandit-Taskar, Steven M. Larson, Nai-Kong V. Cheung. doi:10.1007/s00259-020-05050-z](https://doi.org/10.1007/s00259-020-05050-z)  
* Specifically, the supplemental material of the manuscript contains the technical details for the model and code.
* Please reach out if you have questions or would like to learn more.  
* The modeling project spanned Jul 2017 - Sep 2020.  
* All code is in MATLAB®.


## Notes, Warnings, and Potential Sources of Confusion  
* Omburtamab was originally called 8H9. 8H9 is the name more commonly used in the code.  
* 8H9 is the engineered antibody while B7-H3 is the corresponding natural antigen on tumor cells. (Likewise for antibody 3F8 and antigen GD2).  
* Therapeutic index (TI) is synonymous with therapeutic ratio (TR), which is an older term, but more frequently used in the code.  
* tumorload and R0 both represent the per-volume concentration of antigens. They only differ in their units. tumorload is mol/L, while R0 is antigens/mL. Antigen-density, NR, is related, but is in antigens/cell, which is the raw unit measured. Through calculations and estimations of cell size, you get tumorload and R0 from NR.  
* If figure files are not be created when you run the code as-is, ensure that the appropriate folder structure is in place. The MATLAB® code will create the figure files, but will not create folders if necessary to get to the path that the figure files are supposed to be in.  
* Generally, SI units are used throughout the code, with the exception of mCi (millicurie) being often used as the unit of radioactivity. However, the output of the code (i.e. figures) use MBq (megabecquerel) as the unit of radioactivity.


## Instructions for Understanding the Model's MATLAB® Code  
* The differential equations themselves are provided in _ModOde.m_. 
* However, I recommend running _solveModel.m_ instead, which is a wrapper function that does a multitude of extra pre-processing and post-processing steps (for example calculating AUC). This file is located in the +Simulations folder.
* _DefaultValues.m_ also in the +Simulations folder is meant to be run before any other code is run because it defines the default variables (which can be customized in other files, but should not be changed in _DefaultValues.m_). 
* _Dosage.mlx_ file is a .mlx file I used to create the dosage response figures and are provided as a good example to learn how to run the _solveModel.m_ and _DefaultValues.m files_.  
* Similar to _solveModel.m_, _solveModelMultInfus.m_ is another wrapper function for _ModOde.m_, except that it allows input of multiple infusion times with varying infusion amounts. _SplitDosage.mlx_ file is a .mlx file I used to create the dose fractionation figures and are provided as a good example to learn how to run the _solveModelMultInfus.m_ (it is also another chance to see how _DefaultValues.m_) is used).  

 

## Tested OS and History Details  
* Language used is MATLAB® Version 2017a on a personal Mac computer (macOS Mojave, 2018 MacBook Pro 15in).  
* Originally, I used a Windows computer when working at the MSKCC building in NYC, and then an earlier personal MacBook Pro, before this 2018 personal Mac computer. MATLAB® version might have also been different earlier.  
