# README
![Figure](schematic.png)
[https://excalidraw.com/#json=XBnibVJopP0Gr8kOZ-69F,DMmLNTof7_3ff3TVo4K_0Q]
<!-- [Figure](https://excalidraw.com/#json=XBnibVJopP0Gr8kOZ-69F,DMmLNTof7_3ff3TVo4K_0Q) -->

## The Paper
<!-- Entanglement-Witness-by-Quantum -->

[`Towards efficient and generic entanglement detection by machine learning`](https://arxiv.org/pdf/2211.05592.pdf).
Jue Xu and Qi Zhao, 2022

[![https://arxiv.org/abs/2211.05592](https://img.shields.io/badge/paper%20%28v1%29-arXiv%3A2211.05592-B31B1B)](https://arxiv.org/abs/2211.05592)
<!-- Copyright by Jue Xu (juexu@cs.umd.edu) -->
<!-- [arXiv: 2211.05592](http://arxiv.org/abs/2211.05592). -->
<!-- Entanglement Detection Beyond Measuring Fidelities
M. Weilenmann, B. Dive, D. Trillo, E. A. Aguilar and M. Navascués 
20th December 2019Entanglement Detection Beyond Measuring Fidelities -->


<!-- ####################### INSTALLATION ########################## -->
<!-- ## Installation

The contents of the package do not need any installation beyond the packages required below. -->
<!-- The package can simply be downloaded and run directly.  -->
<!-- Note that the `quantum_state_utils.py` module must be in the same directory as the notebook `Entangle.ipynb` when they are run. -->

## Packages required for numeric simulation: 
- anaconda
- python, version > 3.8.8
- numpy, version > 1.21.5
- matplotlib, version > 3.5.1 
- QuTiP, version > 4.7.0
- Scikit-learn, version > 1.1.2

### [QuTiP](https://qutip.org/)
Installation instructions available at (https://qutip.org/docs/latest/installation).
Reference:
- [basic](https://qutip.org/docs/latest/guide/guide-basics.html)
- [measurement](https://qutip.org/docs/latest/apidoc/functions.html#measurement), https://qutip.org/docs/latest/guide/guide-measurement.html
- [random states](https://qutip.org/docs/latest/apidoc/functions.html#module-qutip.random_objects)
- gates and circuits: https://qutip.org/docs/latest/apidoc/functions.html?#quantum-information-processing, https://qutip.org/docs/latest/guide/qip/qip-simulator.html?#operator-level-circuit-simulation, https://qutip-qip.readthedocs.io/en/stable/qip-basics.html
- [visualization of states](https://qutip.org/docs/latest/apidoc/functions.html?#qutip.visualization.matrix_histogram)

### [Scikit-learn](https://scikit-learn.org/stable/)
Installation instructions available at (https://scikit-learn.org/stable/install).
Reference:

- [SVM](https://scikit-learn.org/stable/modules/svm.html?#support-vector-machines)
and [Kernel for SVM](https://scikit-learn.org/stable/modules/svm.html?#kernel-functions)
- feature elimination: [sklearn.feature_selection.RFE](https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html), [feature_selection](https://scikit-learn.org/stable/auto_examples/feature_selection/plot_rfe_digits.html#sphx-glr-auto-examples-feature-selection-plot-rfe-digits-py)

<!-- #### Recursive feature elimination 
https://stackoverflow.com/questions/71040368/feature-selection-not-working-in-svr-with-rbf-kernel-for-n-features-to-select
https://towardsdatascience.com/powerful-feature-selection-with-recursive-feature-elimination-rfe-of-sklearn-23efb2cdb54e
https://stackoverflow.com/questions/41592661/determining-the-most-contributing-features-for-svm-classifier-in-sklearn -->

## Related projects with open-source code
Classical shadow: [hsinyuan-huang/predicting-quantum-properties](https://github.com/hsinyuan-huang/predicting-quantum-properties); [charleshadfield/adaptiveshadows](https://github.com/charleshadfield/adaptiveshadows)

Unfaithful states: [BenDive/Entanglement-Detection-Beyond-Measuring-Fidelities](https://github.com/BenDive/Entanglement-Detection-Beyond-Measuring-Fidelities)

Neural Network for entanglement classification: [Matlab Code](https://figshare.com/articles/dataset/Codes_and_data_set_of_Transforming_Bell_s_Inequalities_into_State_Classifiers_with_Machine_Learning_/6231662)

<!-- ######################### USAGE ############################# -->

## Usage of code

The [Code](https://github.com/Jue-Xu/Entanglement-Witness/tree/main/Code) folder includes the code supporting the [paper](https://github.com/Jue-Xu/Entanglement-Witness/blob/main/ew.pdf).
The script and notebook are organized as follows. 
[quantum_state_utils.py](https://github.com/Jue-Xu/Entanglement-Witness/blob/main/Code/quantum_state_utils.py) contains the functions that generate specific data (e.g. bi-separable, genuine entangled quantum states) and construct the training/testing datasets. 
The `.py` script must be in the same directory as the notebook when they are run.
The [entangle_detection.ipynb](https://github.com/Jue-Xu/Entanglement-Witness/blob/main/Code/entangle_detection.ipynb) notebook contains the scripts that generates a particular example or a figure.
The [classical_shadow.ipynb](https://github.com/Jue-Xu/Entanglement-Witness/blob/main/Code/classical_shadow.ipynb) notebook demonstrates the performance of classical shadow methods.
<!-- Each of the other files have the scripts that correspond to
a particular example, figure, or table in the manuscript, and can be run independently. -->


<!-- ########################## LICENSE ############################ -->
## License

The code here is made available under the [Creative Commons Attribution license (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)