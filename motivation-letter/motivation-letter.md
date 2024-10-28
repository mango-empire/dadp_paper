---
output: pdf_document
fontsize: 12pt
---

\thispagestyle{empty}
\today

Editor   
The R Journal  
\bigskip

Dear Dr. Mark van der Loo,
\bigskip

Please consider the enclosed manuscript for publication in the R Journal:

Title: "dapper: Data Augmentation for Private Posterior Estimation in R"   
Authors: Kevin Eng, Jordan Awan, Nianqiao Phyllis Ju, Vinayak A. Rao, and Ruobin Gong

Recent advances in data privacy protection such as differential privacy can involve a combination
of deterministic and random transformation to the original data set. These transformations
often invalidate standard statistical methods which typically assume the input data is not
randomly perturbed. While theoretical development of privacy methods have 
moved quickly, their adoption as been relatively slower in part because of a lack of software tools
for analyzing privatized data.

The enclosed manuscript introduces the new package called `dapper` which we hope
will facilitate the adoption of modern differential privacy methods by providing users 
with a flexible tool to perform valid Bayesian inference on data protected by 
differential privacy, allowing them to properly account for the noise 
introduced for privacy protection in their statistical analysis. To the best of our knowledge,
this is the first CRAN R package to provides general purpose tools for conducting Bayesian inference
with privatized data.

We believe the readers of the R Journal will benefit from this article because they
may not be familiar with working with privatized data or if they are, the Bayesian data augmentation
framework used in the manuscript. In particular, the latter requires non-trivial modifications to an existing Bayesian workflow. 
The manuscript goes over key technical aspects that we think
a user would need to know in order to successfully adopt `dapper` to their problem.
The included examples provide a thorough discussion on the model setup which should
serve as a solid template for users to adopt to their own needs.

\bigskip
\bigskip

Regards,
    
Kevin Eng  
Department of Statistics   
Rutgers University  
Piscataway, New Jersey  
ke157@stat.rutgers.edu
