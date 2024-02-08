# model-crossings
Code to run simulations and replicate figures reported in Revsine, Gonzalez-Castillo, Merriam, Bandettini, and Ramirez (2023)

You can find our preprint here: https://www.biorxiv.org/content/10.1101/2023.02.08.527219v1.full

GETTING STARTED: Detailed instructions to run simulatons are included in the m-file named mainMakeFigs.m. We suggest you get started by opening this m-file and readng the included instructions, the summary and warning below, and finally referring to our paper for further details.

# Summary 
This repository consists of matlab functions to achieve two main goals: 

GOAL 1: Conduct simple analyses of low-level image properties of full-images as well as half-images of our two target databases: KDEF and RaFD. Please, request permission to downloand these copyrighted databases here https://rafd.socsci.ru.nl/RaFD2/RaFD?p=main and here https://kdef.se

We focus on the analysis of images of faces with neutral expressions, each photographed in the same five viewpoints (-90 deg to +90 deg, step 45 deg). 

GOAL 2: First generate randomly-connected multi-layered networks considering two network hemispheres and next pass the images of the two relevant databases through the model. Importantly, the instantiated architecture also considers magnification of the central portion of the visual field (analogous to cortical magnification of the fovea in primates) as well as increasing levels of connectivity of units across the two network hemispehres as a function of their hierarchical level along the multi-layered network (emulating the hierarchical arragement of visual areas along the primate ventral visual processing stream). Various pattern analyses are next conducted on the activation patterns associated with the input images in each network layer. See mainMakeFigs.m for details.

# Warning 
Before running mainMakeFigs.m you need to request permission to the copyright holders of the two image databases (web-links provided above). Next, after you download the KDEF and Radboud (RaFD) face databases, please place the associated directories in the locations specified in mainMakeFigs.m and proceed to define the necessary matlab paths (as also explained in mainMakeFigs.m). If the directory structure agrees with the instructions, and the paths to the images and code are correctly defined, running mainMakeFigs.m will replicate all simulations and low-level image feature analyses reported in the paper, as well as reproduce the panels of their associated figures in the paper. Please note that not all figures generated here are shown in the paper. In particular, only a subset of the figures generated in the section of the main-file "Make Figure 7" is shown in the paper. The relevant sub-set of figures can be readily identified by reading the title of each figure.   

# Note 
One of the model variants supported relies on the output of the S1-level representation of the HMAX model as described in Serre et al (2004) and implemented in Matlab by Stanley Bileschi based on the C1 code by Max Riesenhuber and Thomas Serre.


