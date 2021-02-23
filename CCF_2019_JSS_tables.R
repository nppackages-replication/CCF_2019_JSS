########################################################################
## REPLICATION FILES: Calonico, Cattaneo and Farrell (2018)
## nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference
## Last update: 29-April-2019
########################################################################
## This script replicates the simulation tables for different choices of p and deriv

rm(list=ls(all=TRUE))
library(Rcpp)
library(Hmisc)
library(locfit)
library(nprobust)
library(locpol)
library(lpridge)
library(np)

setwd("C:/Users/nsc19/Dropbox/00000_shared--Calonico-Cattaneo-Farrell_2017_nprobust/simuls/CCF_2018_replication")
tryCatch(dir.create("output"))

p     <- 1
deriv <- 0

source("CCF_2019_JSS_simuls.R")
   

