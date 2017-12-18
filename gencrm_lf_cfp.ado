*! v2.0.0, S Bauldry, 18dec2017

capture program drop gencrm_lf_cfp
program gencrm_lf_cfp
  version 15
  
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf xb_c xb_p"
  
  forval i = 1/$nCatm1 {
    local arguments "`arguments' xb_f`i'"
  }
  
  forval i = 1/$nCatm1 {
    local arguments "`arguments' t`i'"
  }
  
  forval i = 2/$nCatm1 {
	local arguments "`arguments' f`i'"
  }
   
  args `arguments'
	
  * tempvars for thresholds and factors
  forval j = 1/$nCatm1 {
    tempvar tau`j'
    qui gen double `tau`j'' = `t`j''
  }
  
  forval j = 2/$nCatm1 {
	tempvar phi`j'
    qui gen double `phi`j'' = `f`j''
  }
  
  * setting values for y
  forval j = 1/$nCat {
    local y_`j' ${y_`j'}
  }
  local M $nCat
	
  *** likelihood function for logit link
  if ( "$Link" == "logit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(invlogit(`tau1' - `xb_c' - `xb_f1' - `xb_p')) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    if ( $nCat == 3 ) {
	  qui replace `lnf' = ln(1 - invlogit(`tau1' - `xb_c' - `xb_f1' - `xb_p')) +  ///
		                  ln(    invlogit(`tau2' - `xb_c' - `xb_f2' - `xb_p'*`phi2')) if $ML_y == `y_2'
	}
	
	if ( $nCat > 3 ) {
	  forval k = 2/$nCatm1 {
	    local meqn_a `" ln(1 - invlogit(`tau1' - `xb_c' - `xb_f1' - `xb_p')) + "'
        local meqn_c `" ln(    invlogit(`tau`k'' - `xb_c' - `xb_f`k'' - `xb_p'*`phi`k'')) "'
    
	    local meqn_b ""
	    local m = `k' - 1
        forval n = 2/`m' {
          local meqn_b `" `meqn_b' ln(1 - invlogit(`tau`n'' - `xb_c' - `xb_f`n'' - `xb_p'*`phi`n'')) + "'
        }
	
        local meqn `" `meqn_a' `meqn_b' `meqn_c' "'
        qui replace `lnf' = `meqn' if $ML_y == `y_`k''
      }
	}
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(`tau1' - `xb_c' - `xb_f1' - `xb_p')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(`tau`o'' - `xb_c' - `xb_f`o'' - `xb_p'*`phi`o'')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }
  
end


/* History
1.0.0  11.22.16  initial likelihood program for arbitrary number of categories
1.1.0  08.25.17  generalized program for unlimited number of categories
1.2.0  09.18.17  fixed bug with non-standard values for Y
2.0.0  12.18.17  new program name, fixed problem with thresholds
