Thank you for these and earlier comments and suggestions. I have expanded the text (see below for point-by-point responses) and also used this as an opportunity submit an updated package version to CRAN (bug fixes, additional utility functions, and updates vignettes based on commens and suggestions).

**line 16/17 how the approach violates the iid assumption is not clear to me, I would appreciate more explanation on this part.**

Good point! I have expanded to explain how the order and prior perception violate the assumption:
"Typically, time series for such multistable stimuli are fitted using Gamma distribution (Fig. 1C). This assumes that individual dominance phase durations are exchangeable, i.e., the order in which they are drawn is unimportant. However, this ignores the effect prior perception through adaptation [@VanEe2009] that makes individual dominance phases serially dependent on prior ones, violating assumptions about independent and identically distributed samples. In other words, each individual dominance phase is drawn from its specific distribution whose parameters are also determined by the prior perceptual history."


**line 18/19: although it cites the prior work on the homogenous fist order process, I would believe that including an explicit model formulation should be easier for readers to see roughly how the underlying model works. Especially in your case, you used Stan, a full Bayesian procedure, so including a simple description for your likelihood function and prior specification would be helpful for users. Or maybe I could suggest linking your vignette page (e.g., https://cran.rstudio.com/web/packages/bistablehistory/vignettes/cumulative-history.html). btw, although it is not needed for the sake of publishing this paper, a better documentation for the model formulation would be more useful (e.g.,image in your vignette, $A$ and $t^prime$, for example, are not defined.**

Yes, thank you! I have included the new section on that. effectively, a large chunk of the corresponding vignette.

**line 26/27: "thus it provides .... ability to compare models via information criteria, etc", this sounds vague to me as frequentist likelihood inference can also do model comparisons using information criteria (e.g., AIC) so maybe better to say explicitly you are referring here WAIC and/or DIC.**

Yes, model comparison via information criteria is not uniquely Bayesian. I have expanded on both the technical details and added a reference on why LOOCV and WAIC are preferred for multilevel models.

**Is there any current R packages or any other software packages for history-dependent analysis? if so, it would be good to cite.**

I have included the information on general time-series packages into the Statement of need.
