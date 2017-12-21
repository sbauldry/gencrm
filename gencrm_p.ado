*! v1.0.0, S Bauldry, 21dec2017

*** This is a placeholder for a predict command and overrides the mistaken
*** predict that comes from ml_p

capture program drop gencrm_p
program define gencrm_p
  dis as error "predict is not currently supported for gencrm"
  exit
end


/* History
1.0.0  12.21.17  placeholder command to override incorrect ml_p
