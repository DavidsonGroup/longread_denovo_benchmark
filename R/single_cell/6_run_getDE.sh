#!/bin/bash
module load R/4.2.3

# DGE (first monocyte, secodn lymphocyte)
Rscript --vanilla dge_limmavoom.R bambu_bambu_oarfish 0,1,6,7,9,10,11 2,3,4,5,8

Rscript --vanilla dge_limmavoom.R isonform_corset_oarfish 0,1,6,8,9,10,11 2,3,4,5,7,12
Rscript --vanilla dge_limmavoom.R isonform_isonclust_oarfish 0,1,7,8,9,11,12 2,3,4,5,6,10

Rscript --vanilla dge_limmavoom.R rattle_corset_oarfish 0,1,6,8,9,10,11 2,3,4,5,7
Rscript --vanilla dge_limmavoom.R rattle_rattle_oarfish 1,2,4,8,9,10,11,12 0,3,5,6,7

Rscript --vanilla dge_limmavoom.R rnabloom2_corset_oarfish 0,1,8,9,10,11,12 2,3,4,5,6,7,13

# DTE DTU
Rscript --vanilla dtu_limmavoom.R bambu_bambu_oarfish 0,1,6,7,9,10,11 2,3,4,5,8

Rscript --vanilla dtu_limmavoom.R isonform_corset_oarfish 0,1,6,8,9,10,11 2,3,4,5,7,12
Rscript --vanilla dtu_limmavoom.R isonform_isonclust_oarfish 0,1,7,8,9,11,12 2,3,4,5,6,10

Rscript --vanilla dtu_limmavoom.R rattle_corset_oarfish 0,1,6,8,9,10,11 2,3,4,5,7
Rscript --vanilla dtu_limmavoom.R rattle_rattle_oarfish 1,2,4,8,9,10,11,12 0,3,5,6,7

Rscript --vanilla dtu_limmavoom.R rnabloom2_corset_oarfish 0,1,8,9,10,11,12 2,3,4,5,6,7,13


