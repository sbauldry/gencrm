*! v2.0.1, S Bauldry, 21dec2017

capture program drop gencrm_lf_n
program gencrm_lf_n
  version 15
		
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf"
		
  forval i = 1/$nCatm1 {
    local arguments "`arguments' xb_f`i'"
  }
  
  args `arguments'
  
  * setting values for y
  forval j = 1/$nCat {
    local y_`j' ${y_`j'}
  }
  local M $nCat
  	
  *** likelihood function for logit link
  if ( "$Link" == "logit" ) {	
    
	* equation for first value of Y
	qui replace `lnf' = ln(invlogit(`xb_f1')) if $ML_y1 == `y_1'
		
	* build equations for middle value of Y
	forval k = 2/$nCatm1 {
      local meqn_b `" ln(invlogit(`xb_f`k'')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - invlogit(`xb_f`n'')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y1 == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(`xb_f1')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(`xb_f`o'')) "'
	}
	qui replace `lnf' = `eqn' if $ML_y1 == `y_`M''
  }



  *** likelihood function for probit link
  if ( "$Link" == "probit" ) {	
    
	* equation for first value of Y
	qui replace `lnf' = ln(normal(`xb_f1')) if $ML_y1 == `y_1'
		
	* build equations for middle value of Y
	forval k = 2/$nCatm1 {
      local meqn_b `" ln(normal(`xb_f`k'')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - normal(`xb_f`n'')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y1 == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - normal(`xb_f1')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - normal(`xb_f`o'')) "'
	}
	qui replace `lnf' = `eqn' if $ML_y1 == `y_`M''
  }



  *** likelihood function for cloglog link
  if ( "$Link" == "cloglog" ) {	
    
	* equation for first value of Y
	qui replace `lnf' = ln(1 - exp(-exp(`xb_f1'))) if $ML_y1 == `y_1'
		
	* build equations for middle value of Y
	forval k = 2/$nCatm1 {
      local meqn_b `" ln(1 - exp(-exp(`xb_f`k''))) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(exp(-exp(`xb_f`n''))) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y1 == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(exp(-exp(`xb_f1'))) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(exp(-exp(`xb_f`o''))) "'
	}
	qui replace `lnf' = `eqn' if $ML_y1 == `y_`M''
  }
  
end


/* History
2.0.0  12.18.17  new program for baseline model
2.0.1  12.21.17  added probit and cloglog

