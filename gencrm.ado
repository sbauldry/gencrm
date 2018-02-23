*! v2.0.7, S Bauldry, 30jan2018

capture program drop gencrm
program gencrm, properties(swml svyb svyj svyr mi or hr eform)
	version 14
	if replay() {
		if (`"`e(cmd)'"' != "gencrm") error 301
		Replay `0'
	}
	else Estimate `0'
end


capture program drop Estimate
program Estimate, eclass sortpreserve
	syntax varlist(numeric fv) [if] [in] [pweight fweight aweight iweight] [, /// 	       
		   FActor(varlist fv)  /// vars with proportionality constraint
		   FRee(varlist fv)    /// vars with no constraint
		   LINK(string)        /// link function (logit, probit, or cloglog)
		   Level(cilevel)      /// display option for confidence intervals
		   vce(string)         /// robust and cluster robust standard errors
		   or hr EForm         /// odds ratio, hazard ratio, and exponential form
		   noLOg               /// -ml model- options
		   svy *               /// -mlopts-, display options
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
	
	* exponential form
	local eform `or' `hr' `eform'
	local efopt : word count `eform'
	if (`efopt' > 1) {
		dis as error "only one of or, hr, eform can be specified"
		exit 198
	}
	
	* check for link function (default to logit)
	if ( "`link'" == "logit" | "`link'" == "l" | "`link'" == "" ) {
		global Link       "logit"
		local  link_title "Ordered Logit Estimates"
		
		if ("`hr'" != "") {
			dis as error "hr not appropriate for logit model"
			exit 198
		}
	}
	else if ( "`link'" == "probit" | "`link'" == "p" ) {
		global Link       "probit"
		local  link_title "Ordered Probit Estimates"
		
		if (`efopt' > 0) {
			dis as error "eform not appropriate for probit model"
			exit 198
		}
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
	if ( "`factor'" != "" ) {
		fvexpand `factor'
		local prIV `r(varlist)'
		
		* verify subset of IVs
		local prchk : list local(prIV) - local(cnIV)
		if ( "`prchk'" != "" ) {
			dis ""
			dis as error "factor{yellow}(`prchk'){red} is not included in" ///
			             " the list of independent variables: {yellow}`IV'"
			exit 198
		}
		
		* remove from list of IVs
		local cnIV : list local(cnIV) - local(prIV)
	}
	

	* create ML model statements
	
	* case 0: baseline model with no covariates
	if ( "`free'" == "" & "`factor'" == "" & "`cnIV'" == "" ) {
		forval i = 1/$nCatm1 {
			local model "`model' (tau`i': `Y' = )"
		}
		
		* obtain ML estimates
		ml model lf gencrm_lf_n `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing maximize
			
		* replace current b, V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
		local eqno = 0
	}
	
	
	* case 1: all variables with parallel assumption
	if ( "`free'" == "" & "`factor'" == "" & "`cnIV'" != "" ) {
		local model "(parallel: `Y' = `cnIV', nocons)"
	
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
	
		* obtain ML estimates
		local eqno = 1
		ml model lf gencrm_lf_c `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
			
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 2: all variables with non-parallel assumption
	if ( "`free'" != "" & "`factor'" == "" & "`cnIV'" == "" ) {
		local model "(eq1: `Y' = `free', nocons)"
		
		forval i = 2/$nCatm1 {
			local model "`model' (eq`i': `free', nocons)"
		}
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		* obtain ML estimates
		local eqno = $nCatm1
		ml model lf gencrm_lf_f `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
		
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 3: all variables with proportionality assumption
	if ( "`free'" == "" & "`factor'" != "" & "`cnIV'" == "" ) {
		local model "(factor: `Y' = `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}
		
		* obtain ML estimates
		local eqno = 1
		ml model lf gencrm_lf_p `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
		
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 4: subset of variables constrained and proportionality assumption
	if ( "`free'" == "" & "`factor'" != "" & "`cnIV'" != "" )  {
		local model "(parallel: `Y' = `cnIV', nocons) (factor: `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}	
		
		* obtain ML estimates
		local eqno = 2
		ml model lf gencrm_lf_cp `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
		
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 5: subset of variables constrained and non-parallel assumption
	if ( "`free'" != "" & "`factor'" == "" & "`cnIV'" != "" )  {
		local model "(parallel: `Y' = `cnIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' (eq`i': `free', nocons)"
		}
		
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
		
		* obtain ML estimates
		local eqno = $nCat
		ml model lf gencrm_lf_cf `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
		
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 6: subset of variables proportionality and non-parallel assumption
	if ( "`free'" != "" & "`factor'" != "" & "`cnIV'" == "" )  {
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
		local eqno = $nCat
		ml model lf gencrm_lf_fp `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
		
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 7: subset of variables with non-parallel assumption and 
	*         proportionality assumption
	if ( "`free'" != "" & "`factor'" != "" & "`cnIV'" != "" )  {
		local model "(parallel: `Y' = `cnIV', nocons) (factor: `prIV', nocons)"
		
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
		local eqno = $nCat + 1
		ml model lf gencrm_lf_cfp `model' `wgt' if `touse', ///
		  title(`link_title') vce(`vcetype') `log' `mlopts' missing ///
		  waldtest(`eqno') maximize
		
		* replace current b and V
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}	
	
	* return and display results
	ereturn scalar k_cat = $nCat
	ereturn scalar k_eform = `eqno'
	ereturn local cmd gencrm
	ereturn local free `free'
	ereturn local factor `factor'
	ereturn local link `link'
	ereturn local vce "`vce'"
	ereturn local vceptype "`vcetype'"
	ereturn local clustvar "`clustervar'"
	ereturn local predict "gencrm_p"
	ereturn repost b = `b' V = `v', rename 
	
	Replay, level(`level') `eform' `diopts' `options'
end



capture program drop Replay
program Replay
	syntax [, Level(cilevel) or hr EForm *]
	
	* display options
	_get_diopts diopts options, `options'
	local diopts `diopts' `or' `hr' `eform' level(`level') 
	
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
2.0.2  12.20.17  added nolog as default option
2.0.3  12.21.17  updated ML options
2.0.4  12.21.17  removed eform option and set predict
2.0.5  12.21.17  fixed eform options
2.0.6  01.15.18  fixed bug with Wald test
2.0.7  01.30.18  fixed eform options for inappropriate links
2.0.8  02.18.18  changed name of prop option to factor option
2.0.9  02.18.18  set command to work with Stata v14
2.0.10 02.23.18  relabeled constrained as parallel

