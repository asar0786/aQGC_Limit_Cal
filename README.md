# Info

- The two codes `WVChannel_GetCard.C` and `WVChannel_GetCard_WithHiggsDistributions.C` are same except the first one does not contain doubly charged higgs while the second one contains.
- If you get **segmentation fault** with code `WVChannel_GetCard_WithHiggsDistributions.C` then check if the doubly charged higgs sample info are not commented in the sample text file (e.g. this file [DibosonBoostedElMuSamples13TeV_WWTree_CommonNtuple_For1and2Lepton_MuonPtScale_2018_07_24_10h36.txt](https://github.com/ram1123/ControlPlots/blob/master/LimitCardGenerate/DibosonBoostedElMuSamples13TeV_WWTree_CommonNtuple_For1and2Lepton_MuonPtScale_2018_07_24_10h36.txt)).


# Prepare Input RootFile For Combine

> root -l -b -q WVChannel_GetCard.C

> root -l -b -q ZVChannel_GetCard.C

Output root file from background estimation and the associated all plots for alpha and up/dn can be found in root file present here with name: [WV_bkg_estimation_4Bins.root](https://github.com/ram1123/ControlPlots/blob/master/LimitCardGenerate/WV_bkg_estimation_4Bins.root) and [ZV_bkg_estimation_4Bins.root](https://github.com/ram1123/ControlPlots/blob/master/LimitCardGenerate/ZV_bkg_estimation_4Bins.root)
