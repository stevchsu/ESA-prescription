# Code Description
This repository contains the code necessary to generate the main results described in our paper "Computer-assisted prescription of Erythropoiesis-stimulating agents (ESA) in Patients under maintenance hemodialysis, a randomized clinical trial for AI model selection".

The structure of the repository is as follows:

- In `R/prepare_data.R`, we preprocess data to prepare data for both bagging REEM trees and MERF.

- In the folder `R/model`, there are codes of our four models, while LSTM-1 and LSTM-2 share the same network architecture but toward different target Hbc levels. We also infer ESA dose with the two LSTM with `R/model/LSTM_1_and_2.R`.

- In the folder `R/dose_inference`, there are codes to make inferences on ESA doses with bagging REEM trees and MERF.

- In the folder `R/lib,` we developed the necessary functions to build MERF based on the R package `ranger`.