# gencrm
This repository contains a set of commands for fitting generalized continuation ratio models in Stata. The command allows for effects to be freely estimated across thresholds, to have a proportionality constraint across thresholds, or to be constrained to be equal across thresholds. Type "net install gencrm, from(https://github.com/sbauldry/gencrm/raw/master) replace" without the quotes in Stata's command window to install the commands.

#Files
1.  gencrm.ado:        primary command file
2.  gencrm_lf_n.ado:   likelihood function for baseline model with no covariates
3.  gencrm_lf_c.ado:   likelihood function for constrained model (parallel lines)
4.  gencrm_lf_f.ado:   likelihood function for unconstrained model (no parallel lines)
5.  gencrm_lf_p.ado:   likelihood function for proportionality constraint model
6.  gencrm_lf_cf.ado:  likelihood function mix of constrained and free
7.  gencrm_lf_cp.ado:  likelihood function mix of constrained and proportionality
8.  gencrm_lf_fp.ado:  likelihood function mix of free and proportionality
9.  gencrm_lf_cfp.ado: likelihood function for general model
10. gencrm_p.ado:      placeholder command for predict
11. gencrm.pkg:        Stata package file
12. stata.toc:         Stata toc file
