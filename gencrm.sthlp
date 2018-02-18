{smcl}
{* Last revised 02/18/18}{…}
{title:Title}

{p2colset 8 18 24 2}{...}
{p2col:{cmd:gencrm}} Generalized Continuation Ratio Regression Models {p_end}
{p2colreset}{...}

{marker syn}
{title:Syntax}

{p 8 18 2}
{cmd:gencrm} {it:depvar} [{it:indepvars}] 
{ifin} {weight} 
[{cmd:,}
{cmdab:fa:ctor(}{it:varlist}{cmd:)}
{cmdab:fr:ee(}{it:varlist}{cmd:)}
{cmdab:link(}{it:string}{cmd:)}
{cmdab:vce(}{it:vcetype}{cmd:)}
{cmdab:or} {cmdab:ef:orm}
{it:display_options} {it:maximize_options}
{p_end}


{title:Description}

{p 4 4 2} {cmd:gencrm} is a user-written command that fits generalized continuation ratio models (Fullerton and Xu 2016). This class of models permits three types of covariates: (1) covariates whose coefficients are constrained to be equal across cutpoint equations, (2) covariates whose coefficients vary across cutpoint equations by a common factor (listed using the prop option), and (3) covariates whose coefficients freely vary across cutpoint equations (listed using the free option). 

{p 4 4 2} {cmd:gencrm} supports factor variables and the {cmdab:mi estimate} prefix. Note that the command is not currently integrated with Stata's {cmdab:predict} or {cmdab:margins} commands. Because of this, the command also is not integrated with Stata's {cmdab:svy} suite of commands, though it is possible to incorporate weights.


{title:Options}

{p 4 8 2} {cmdab:fa:ctor(}{it:varlist}{cmd:)} specifies which, if any, of the independent variables have coefficients that vary across cutpoint equations by a common factor (i.e., proportionality constraint)

{p 4 8 2} {cmdab:fr:ee(}{it:varlist}{cmd:)} specifies which, if any, of the independent variables have coefficients that freely vary across cutpoint equations.

{p 4 8 2} {cmdab:link(}{it:string}{cmd:)} specifies the link function to be used. The default link function is the logit link. Users may also specify a probit or a cloglog link.

{p 4 8 2} {cmdab:vce(}{it:vcetype}{cmd:)} specifies the type of standard error reported. Users may specify {cmdab:r:obust} or {cmdab:cl:uster} {it:clustvar} standard errors. 

{p 4 8 2} {cmdab:or} report odds ratios

{p 4 8 2} {cmdab:ef:orm} report exponentiated coefficients


{title:Examples}
{p 4 4 2}
The following example draws on a Stata dataset with information about birthweight to illustrate a few variants of generalized continuation ratio models. Note that this is not a particularly good substantive example as birthweight is not an outcome that we typically think of as progressing through a series of stages. Nonetheless, it is useful to illustrate the {cmdab:gencrm} command.

{p 4 4 2}
Example 1: {it}Continuation ratio model with parallel lines assumption for all variables.{sf}

{p 8 12 2}{cmd: webuse lbw}{p_end}
{p 8 12 2}{cmd: recode bwt (0/2499 = 1) (2500/2999 = 2) (3000/3499 = 3) (3500/5000 = 4)}{p_end}
{p 8 12 2}{cmd: gencrm bwt smoke i.race lwt ht ui}{p_end}


{p 4 4 2}
Example 2: {it}Continuation ratio model with no parallel lines assumption for smoking.{sf}

{p 8 12 2}{cmd: gencrm bwt smoke i.race lwt ht ui, free(smoke)}{p_end}


{p 4 4 2}
Example 3: {it}Continuation ratio model with proportionality constraint for smoking.{sf}

{p 8 12 2}{cmd: gencrm bwt smoke i.race lwt ht ui, factor(smoke)}{p_end}




{marker aut}
{title:Authors}

{p 4 4 2}
Shawn Bauldry {break}
Purdue University {break}
Department of Sociology {break}
sbauldry@purdue.edu {break}

{p 4 4 2}
Jun Xu {break}
Ball State University {break}
Department of Sociology {break}
jxu@bsu.edu {break}

{p 4 4 2}
Andrew Fullerton {break}
Oklahoma State University {break}
Department of Sociology {break}
andrew.fullerton@okstate.edu {break}



{title:References}

{p 4 8 2}Fullerton, A and Xu, J. 2016. {it:Ordered Regression Models: Parallel, Partial, and Non-Parallel Alternatives}. New York: CRC Press.
