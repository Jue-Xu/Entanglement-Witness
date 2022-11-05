Introduction:

Here we provide codes in order to help readers reproduce main numerical results in our manuscript "Transforming Bell's Inequalities into State Classifiers with Machine Learning".

--------------------------------------------


Preparation:

install Matlab 2017 + python 3.6 or later
Matlab codes generate data, and python codes train, test these data and plot its result.

install QETLAB (http://www.qetlab.com/Main_Page) in Matlab

install keras 2.1.5 with tensorflow 1.8.0 as backend in python

We recommend to install Anaconda for buliding python enviroment, then just input "pip install tensorflow" and "pip install keras" in command line to install the latest version of these packages.
If you have an Nvidia card, see "https://www.tensorflow.org/install/" for more details to install the gpu version of tensorflow.

All of codes have been tested in our ThinkStation P910 (256 GB ARM) with a GTX Titan X card.

--------------------------------------------


Begin:
Extract all rar files (pre-trained results are provided for illustration in "npz" directory), set up an empty directory "mat" and run the codes below for reproducing the results.

--------------------------------------------

Section
"Optimizing CHSH operator with Machine Learning"

Matlab
GenData.main_test_run()

Python (in ./py)
trainer/TrainContinuousHidden.train()
predictor/PredictContinuous(qnum=2)
ploter/PlotContinuous.plot3d_for_paper(selection=(0, 1)) 
# or set selection=(0, 2)

--------------------------------------------

Section
"Classifying general two-qubit states"

Matlab
GenData.main_general_detection() 
% set istrain=1 for train, 0 for test, both data should be generated once

Python
trainer/TrainGap.train()
predictor/PredictGap()

analyzing_hidden_units.ipynb (first open jupyter notebook)

ploter/PlotGapAndBlackBox.plot_all()
# drag the legends to a proper position

--------------------------------------------

Section
"Machine learning for identifying bound entangled states"

Matlab
GenData.main_UPB()

Python
trainer/TrainBisep.train()
ploter/PlotUPB()

--------------------------------------------

Note that the codes only train one model in "Classifying general two-qubit states" and "Machine learning for identifying bound entangled states" because complete training costs much longer time. Just follow the comment in the codes, reset some hyper parameters to training all models in order to plot.

Original generated "quantum states" (mat files) are too large to be uploaded.



