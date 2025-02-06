Title of the dataset:
Fitting of two-step unfolding model to circular dichroism spectra

Creators:
Celia Fricke
Lars Boyens-Thiele

Contributors:
Celia Fricke
Lars Boyens-Thiele

Related publication:
On the thermal and chemical stability of DNAJB6b and its globular domains

Description:
Fitting of circular dichroism spectra at various guanidine hydrochloride concentrations to a two-step unfolding model.

Keywords:
unfolding midpoints
DNAJB6b
Circular dichroism spectroscopy

This dataset contains the following files:
denaturation_spectra_data_0to5M.txt
two_step_unfolding_model.ipynb

Explanation of variables:
s1: cd spectrum of native form
s2: cd spectrum of intermediate form
s3: cd spectrum of fully denatured form
p: array of fitting parameters:
m_I: m-value for first transition
m_D: m-value for second transition
CM_I: Midpoint denaturant concentration of first transition
CM_D: Midpoint denaturant concentration of second transition
conc: Denaturant concentration

Methods, materials and software:
A  three-state unfolding model was fitted to the measured CD spectra, by simultaneously solving for the optimal combination of the two denaturation midpoints, the two m-values, and the spectra of the three pure states. Each measured spectrum was modelled as a linear combination of the spectra of the three states, according to the expected populations from the unfolding model. In addition to minimizing the residual sum of squares, a first-order Tikhonov regularization term was also applied to the three spectra to ensure smooth solutions. Minimization of the cost function was done using the Broyden–Fletcher–Goldfarb–Shanno algorithm as implemented in the minimize function of scipy.optimize (10.1038/s41592-019-0686-2). The spectra measured at 0 M, 2 M, and 5.2 M GdnHCl were used as initial guesses for the three pure spectra.
The fractions of the three species were calculated assuming a linear dependence of the free energy differences on the GdnHCl concentration.

This dataset is published under the BSD 3-Clause license.