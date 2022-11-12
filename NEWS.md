
# Version 0.1.1 (current version on github)

- **(R) formatQM:** Fixed missing `var` argument in formatQM function and in each method when formatQM is called.

# Version 0.1.0

- **General:** Incorporated the capacity to replace low precipitation values (values below `pp_threshold`) with random values between 0 and `pp_factor`. This capacity also solves the problem when the moving window has only cero values (dry months or arid zones), which causes numerical problems with the invers Cumulative Distribution Functions. This can be disabled by defining `pp_threshold = 0`. `pp_threshold` and `pp_factor` are optional arguments defined with values of `1` and `1/100`, respectively.
* **QDM:** Incorporated the capacity to truncate the scaling factor for precipitations in the QDM method to avoid overestimations in dry periods caused by numerical indefinitions of the scaling factor. The scaling factor is truncated to a value defined with the parameter `rel_change_th` when the denominator of the scaling factor is below the parameter `inv_mod_th`. This can be disabled by defining `inv_mod_th = 0`. `rel_change_th` and `inv_mod_th` are optional arguments defined with values of `2` and `pp_threshold`, respectively.
* **(R) Report function:** Fixed the xlims and ylims of the CDF and monthly series plots, respectively.

# Version 0.0.0

First release of the climQMBC package.
