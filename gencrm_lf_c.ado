*! v2.0.1, S Bauldry, 21dec2017

capture program drop gencrm_lf_c
program gencrm_lf_c
  version 15
  
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf xb_c"
  
  forval i = 1/$nCatm1 {
    local arguments "`arguments' t`i'"
  }
	
  args `arguments'
	
  * tempvars for thresholds
  forval j = 1/$nCatm1 {
    tempvar tau`j'
    qui gen double `tau`j'' = `t`j''
  }
  
  * setting values for y
  forval j = 1/$nCat {
    local y_`j' ${y_`j'}
  }
  local M $nCat
	
  *** likelihood function for logit link
  if ( "$Link" == "logit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(invlogit(`tau1' - `xb_c')) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(invlogit(`tau`k'' - `xb_c')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - invlogit(`tau`n'' - `xb_c')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(`tau1' - `xb_c')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(`tau`o'' - `xb_c')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }

  

  *** likelihood function for probit link
  if ( "$Link" == "probit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(normal(`tau1' - `xb_c')) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(normal(`tau`k'' - `xb_c')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - normal(`tau`n'' - `xb_c')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - normal(`tau1' - `xb_c')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - normal(`tau`o'' - `xb_c')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }
  

  
  *** likelihood function for clogog link
  if ( "$Link" == "cloglog" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(1 - exp(-exp(`tau1' - `xb_c’))) if $ML_y == `y_1'
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(1 - exp(-exp(`tau`k'' - `xb_c’))) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(exp(-exp(`tau`n'' - `xb_c’))) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(exp(-exp(`tau1' - `xb_c’))) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(exp(-exp(`tau`o'' - `xb_c’))) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }


end


/* History
1.0.0  11.15.16  initial likelihood program for arbitrary number of categories
1.1.0  08.25.17  generalized program for unlimited number of categories
1.2.0  09.18.17  fixed bug with non-standard values for Y
2.0.0  12.18.17  new program name and fixed likelihood functions
2.0.1  12.21.17  added probit and cloglog
