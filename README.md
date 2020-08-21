# JointMM
A joint modeling framework for analyzing zero-inflated longitudinal proportions and time to event data

This package is developed to investigate the association between zero-inflated longitudinal proportions and time to an event,e.g., disease onset. It is in particular applicable to **microbiome proportion data (or called relative abundance data)** from prospective studies. The model JointMM is specifically designed to handle the zero-inflated and highly skewed longitudinal microbial proportion data to examine whether the temporal pattern of microbial presence and/or the non-zero microbial proportions are associated with differences in the time to an event. 

Installation of JointMM in R:

library("devtools");
install_github("JiyuanHu/JointMM");

Reference: Hu J, Wang C, Blaser M, Li H (2020). Joint modeling of zero-inflated longitudinal proportions and time-to-event data with application to a gut microbiome study. (Submitted)
