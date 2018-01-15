*! v2.0.2, S Bauldry, 15jan2018

capture program drop gencrm_lf_p
program gencrm_lf_p
  version 15
  
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf xb_p"
  
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
    qui replace `lnf' = ln(invlogit(`tau1' - `xb_p')) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    if ( $nCat == 3 ) {
	  qui replace `lnf' = ln(1 - invlogit(`tau1' - `xb_p')) +  ///
		                  ln(    invlogit(`tau2' - `xb_p'*`phi2')) if $ML_y == `y_2'
	}
	
	if ( $nCat > 3 ) {
	  forval k = 2/$nCatm1 {
	    local meqn_a `" ln(1 - invlogit(`tau1' - `xb_p')) + "'
        local meqn_c `" ln(    invlogit(`tau`k'' - `xb_p'*`phi`k'')) "'
    
	    local meqn_b ""
	    local m = `k' - 1
        forval n = 2/`m' {
          local meqn_b `" `meqn_b' ln(1 - invlogit(`tau`n'' - `xb_p'*`phi`n'')) + "'
        }
	
        local meqn `" `meqn_a' `meqn_b' `meqn_c' "'
        qui replace `lnf' = `meqn' if $ML_y == `y_`k''
      }
	}
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(`tau1' - `xb_p')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(`tau`o'' - `xb_p'*`phi`o'')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }



  *** likelihood function for probit link
  if ( "$Link" == "probit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(normal(`tau1' - `xb_p')) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    if ( $nCat == 3 ) {
	  qui replace `lnf' = ln(1 - normal(`tau1' - `xb_p')) +  ///
		                  ln(    normal(`tau2' - `xb_p'*`phi2')) if $ML_y == `y_2'
	}
	
	if ( $nCat > 3 ) {
	  forval k = 2/$nCatm1 {
	    local meqn_a `" ln(1 - normal(`tau1' - `xb_p')) + "'
        local meqn_c `" ln(    normal(`tau`k'' - `xb_p'*`phi`k'')) "'
    
	    local meqn_b ""
	    local m = `k' - 1
        forval n = 2/`m' {
          local meqn_b `" `meqn_b' ln(1 - normal(`tau`n'' - `xb_p'*`phi`n'')) + "'
        }
	
        local meqn `" `meqn_a' `meqn_b' `meqn_c' "'
        qui replace `lnf' = `meqn' if $ML_y == `y_`k''
      }
	}
	
	* build equation for last value of Y
	local eqn `" ln(1 - normal(`tau1' - `xb_p')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - normal(`tau`o'' - `xb_p'*`phi`o'')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }



  *** likelihood function for cloglog link
  if ( "$Link" == "cloglog" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(1 - exp(-exp(`tau1' - `xb_p'))) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    if ( $nCat == 3 ) {
	  qui replace `lnf' = ln(    exp(-exp(`tau1' - `xb_p'))) +  ///
		                  ln(1 - exp(-exp(`tau2' - `xb_p'*`phi2'))) if $ML_y == `y_2'
	}
	
	if ( $nCat > 3 ) {
	  forval k = 2/$nCatm1 {
	    local meqn_a `" ln(    exp(-exp(`tau1' - `xb_p'))) + "'
        local meqn_c `" ln(1 - exp(-exp(`tau`k'' - `xb_p'*`phi`k''))) "'
    
	    local meqn_b ""
	    local m = `k' - 1
        forval n = 2/`m' {
          local meqn_b `" `meqn_b' ln(exp(-exp(`tau`n'' - `xb_p'*`phi`n''))) + "'
        }
	
        local meqn `" `meqn_a' `meqn_b' `meqn_c' "'
        qui replace `lnf' = `meqn' if $ML_y == `y_`k''
      }
	}
	
	* build equation for last value of Y
	local eqn `" ln(exp(-exp(`tau1' - `xb_p'))) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(exp(-exp(`tau`o'' - `xb_p'*`phi`o''))) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }
  
end


/* History
2.0.0  12.18.17  new program for all proportionality constraint
2.0.1  12.21.17  added probit and cloglog
2.0.2  01.15.18  fixed bug with parentheses in cloglog link
