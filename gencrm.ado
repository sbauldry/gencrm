*! v2.0.1, S Bauldry, 19dec2017

capture program drop gencrm
program gencrm, properties(swml svyb svyj svyr mi or eform)
	version 15
	if replay() {
		if (`"`e(cmd)'"' != "gencrm") error 301
		Replay `0'
	}
	else Estimate `0'
end


capture program drop Estimate
program Estimate, eclass sortpreserve
	syntax varlist(numeric fv) [if] [in] [pweight fweight aweight iweight] [, /// 	       
		   PRop(varlist fv)    /// vars with proportionality constraint
		   FRee(varlist fv)    /// vars with no constraint
		   LINK(string)        /// link function (logit, probit, or cloglog)
		   Level(cilevel)      /// display option for confidence intervals
		   vce(string)         /// robust and cluster robust standard errors
		   or EForm            /// exponential form options
		   svy *               /// -mlopts, display options
		   ]
		   
	* sample selection
	marksample touse
		   
	* parse VCE statement
	if ( `"`vce'"' != "" ) {
		my_vce_parse , vce(`vce')
		local vcetype    "robust"
		local clustervar "`r(clustervar)'"
		if "`clustervar'" != "" {
			markout `touse' `clustervar'
		}
	}
	
	* parse weight statement
	if ( "`weight'" != "" ) {
		local wgt "[`weight'`exp']"
		markout `touse' `wvar'
	}
	
	* check that there are cases
	qui count if `touse' != 0
	if r(N) == 0 {
		dis as error "There are no observations"
		exit 2000
	}
	
	* check for link function (default to logit)
	if ( "`link'" == "logit" | "`link'" == "l" | "`link'" == "" ) {
		global Link       "logit"
		local  link_title "Ordered Logit Estimates"
	}
	else if ( "`link'" == "probit" | "`link'" == "p" ) {
		global Link       "probit"
		local  link_title "Ordered Probit Estimates"
	} 
	else if ( "`link'" == "cloglog" | "`link'" == "c" ) {
		global Link       "cloglog"
		local  link_title "Ordered Complementary Log-Log Estimates"
	} 
	else {
		dis ""
		dis as error "{yellow}`link'{red} is not a supported link function"
		exit 198
	}
	
	* exponential form
	local eform `or' `eform'
	local efopt : word count `eform'
	if `efopt' > 1 {
		dis as error "only one of or or eform can be specified"
		exit 198
	}
	
	* check for display and maximization options
	_get_diopts diopts options, `options'
	mlopts mlopts, `options'
		
	* identify DV and IVs
	gettoken Y X : varlist
	_fv_check_depvar `Y'

	* set globals for # categories
	tempname Yval
	qui tab `Y' if `touse', matrow(`Yval') 
	global nCat   = r(r)
	global nCatm1 = r(r) - 1
	
	macro drop y_*
	forval i = 1/$nCat {
		global y_`i' = `Yval'[`i',1]
	}
		
		* too few categories
		if ( $nCat < 3 ) {
			dis ""
			dis as error "{yellow}`Y'{red} has $nCat categories - a minimum" ///
			             " of 3 is required"
			exit 148
		}
		
	* prepare IVs
	if ( "`X'" != "" ) {
		fvexpand `X'
		local cnIV `r(varlist)'
	}
	
	* prepare no constraint variables
	if ( "`free'" != "" ) {
		fvexpand `free'
		local frIV `r(varlist)'
		
		* verify subset of IVs
		local frchk : list local(frIV) - local(cnIV)
		if ( "`frchk'" != "" ) {
			dis ""
			dis as error "free{yellow}(`frchk'){red} is not included in" ///
			             " the list of independent variables: {yellow}`IV'"
			exit 198
		}
		
		* remove from list of IVs
		local cnIV : list local(cnIV) - local(frIV)
	}
	
	* prepare proportionality constraint variables
	if ( "`prop'" != "" ) {
		fvexpand `prop'
		local prIV `r(varlist)'
		
		* verify subset of IVs
		local prchk : list local(prIV) - local(cnIV)
		if ( "`prchk'" != "" ) {
			dis ""
			dis as error "prop{yellow}(`prchk'){red} is not included in" ///
			             " the list of independent variables: {yellow}`IV'"
			exit 198
		}
		
		* remove from list of IVs
		local cnIV : list local(cnIV) - local(prIV)
	}
	

	* create ML model statements
	
	* case 0: baseline model with no covariates
	if ( "`free'" == "" & "`prop'" == "" & "`cnIV'" == "" ) {
		forval i = 1/$nCatm1 {
			local model "`model' (tau`i': `Y' = )"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_n `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
			
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	
	* case 1: all variables with parallel assumption
	if ( "`free'" == "" & "`prop'" == "" & "`cnIV'" != "" ) {
		local model "(constrained: `Y' = `cnIV', nocons)"
	
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
	
		* obtain ML estimates
		ml model lf gencrm_lf_c `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
			
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 2: all variables with non-parallel assumption
	if ( "`free'" != "" & "`prop'" == "" & "`cnIV'" == "" ) {
		local model "(eq1: `Y' = `free', nocons)"
		
		forval i = 2/$nCatm1 {
			local model "`model' (eq`i': `free', nocons)"
		}
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_f `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 3: all variables with proportionality assumption
	if ( "`free'" == "" & "`prop'" != "" & "`cnIV'" == "" ) {
		local model "(factor: `Y' = `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_p `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 4: subset of variables constrained and proportionality assumption
	if ( "`free'" == "" & "`prop'" != "" & "`cnIV'" != "" )  {
		local model "(constrained: `Y' = `cnIV', nocons) (factor: `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}	
		
		* obtain ML estimates
		ml model lf gencrm_lf_cp `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 5: subset of variables constrained and non-parallel assumption
	if ( "`free'" != "" & "`prop'" == "" & "`cnIV'" != "" )  {
		local model "(constrained: `Y' = `cnIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' (eq`i': `free', nocons)"
		}
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_cf `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 6: subset of variables proportionality and non-parallel assumption
	if ( "`free'" != "" & "`prop'" != "" & "`cnIV'" == "" )  {
		local model "(factor: `Y' = `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' (eq`i': `free', nocons)"
		}
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_fp `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 7: subset of variables with non-parallel assumption and 
	*         proportionality assumption
	if ( "`free'" != "" & "`prop'" != "" & "`cnIV'" != "" )  {
		local model "(constrained: `Y' = `cnIV', nocons) (factor: `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' (eq`i': `free', nocons)"
		}
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_cfp `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}	
	
	* Support for margins
        * note: not implemented in this version
	forval i = 1/$nCat {
		local j = `Yval'[`i',1]
		local mdflt `mdflt' predict(pr outcome(`j'))
	}

	* return and display results
	ereturn scalar k_cat = $nCat
	ereturn local cmd cratreg
	ereturn local free `free'
	ereturn local prop `prop'
	ereturn local link `link'
	ereturn local vce "`vce'"
	ereturn local vceptype "`vcetype'"
	ereturn local clustvar "`clustervar'"
	ereturn repost b = `b' V = `v', rename 
	
	Replay, level(`level') `eform' `diopts' `options'
end



capture program drop Replay
program Replay
	syntax [, Level(cilevel) or EForm *]
	
	* display options
	_get_diopts diopts options, `options'
	local diopts `diopts' `eform' level(`level') `or' 
	
	ml display, `diopts'
end



capture program drop my_vce_parse
program my_vce_parse, rclass
	syntax [, vce(string) ]
	
	local case : word count `vce'
	
	if ( `case' > 2 ) {
		dis `"{red}{bf:vce(`vce')} invalid"'
		exit 498
	}
	
	local 0 `", `vce'"'
	syntax [, Robust Cluster * ]
	
	if ( `case' == 2 ) {
		if "`robust'" == "robust" | "`cluster'" == "" {
			dis `"{red}{bf:vce(`vce')} invalid"'
			exit 498
		}
		
		capture confirm numeric variable `options'
		if _rc {
			dis `"{red}{bf:vce(`vce')} invalid"'
			exit 498
		}
		local clustervar "`options'"
	}
	else {
		if ( "`robust'" == "" ) {
			dis `"{red}{bf:vce(`vce')} invalid"'
			exit 498
		}
	}
	
	return clear
	return local clustervar "`clustervar'"
end


/* History
1.0.0  11.15.16  initial program for arbitrary number of categories
1.0.1  06.21.17  updated labels
1.0.2  07.16.17  updated labels again
1.1.0  08.29.17  changed name of program
1.2.0  09.18.17  fixed bug with non-standard values for Y
1.2.1  09.20.17  updated exponential form options
2.0.0  12.18.17  new program name and fixed likelihood functions
2.0.1  12.19.17  added baseline model with no covariates

