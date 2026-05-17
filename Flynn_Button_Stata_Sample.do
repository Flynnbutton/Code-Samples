********************************************************************************
*                                                                              *
*       Library Access and School Test Scores - KS4 Analysis (Excerpt)        *
*                                                                              *
********************************************************************************
*
* This is an excerpt from a larger do-file developed as part of a Master's
* research project at the Barcelona School of Economics. The full script
* includes treatment variable construction, main regressions and robustness
* checks, event study specifications with pre-trend F-tests, heterogeneity
* analysis by school disadvantage quintile, data visualisation, and
* publication-ready regression tables exported to LaTeX via esttab.
*
* This excerpt covers:
*   Section 1 - Data setup and cleaning
*   Section 2 - Within-year standardisation and measure transition validation
*   Section 3 - Treatment variable construction
*   Section 4 - Additional controls and flags
*   Section 5 - Main regressions
*
* Author: Flynn Button. Developed as part of a collaborative graduate research
* project at the Barcelona School of Economics.
*
********************************************************************************


clear all

global root    "/Users/Flynn/Desktop/BSE/T3/Masters Project/school_panel_masters_v9"
global data    "$root/Data"
global code    "$root/Code"
global outputs "$root/Outputs"


********************************************************************************
**#                         1. Data Setup                                     **
********************************************************************************

import delimited "${data}/school_panel_ks4_v9.csv", clear

sort urn year

* Destring Variables
destring totpups tpup att8scr att8screng att8screbac att8scrmat att8scropen ///
    ks2aps la_labpart la_unemp la_earnings ttaps pt5em_94 ptebaceng_94 ///
    ptebaceng_e_ptq_ee ptebacmat_94 ptebacmat_e_ptq_ee ///
    la_tot_exp la_educ_exp la_cultural_exp, replace force

destring ptfsm6cla1a ptealgrp2, replace ignore("%") force
replace ptealgrp2 = ptealgrp2/100

* Drop Special Schools and Covid Years
drop if inlist(minorgroup, "Special School", "Special school")
drop if year == 2020 | year == 2021


********************************************************************************
**#                   2. Within-Year Standardisation                          **
********************************************************************************

capture drop mean_ttaps sd_ttaps z_ttaps mean_att8scr sd_att8scr z_att8scr


* Standardise both attainment measures within year
foreach var of varlist ttaps att8scr {
    bysort year: egen mean_`var' = mean(`var')
    bysort year: egen sd_`var'  = sd(`var')
    gen z_`var' = (`var' - mean_`var') / sd_`var'
}

* Combined standardised outcome: ttaps pre-2016, att8scr post-2015
* (England changed attainment measure in 2016)
gen overall_std = z_ttaps   if year <= 2015
replace overall_std = z_att8scr if year >  2015

* Additional standardised outcome variables
egen pct5_std   = std(pt5em_94)
egen pcteng_std = std(ptebaceng_94)
egen pctmath_std = std(ptebacmat_94)



* Measure Transition Validation
* Check rank correlation between old and new measures at the transition point
xtile q_old = ttaps   if year == 2015, nquantiles(5)
xtile q_new = att8scr if year == 2016, nquantiles(5)

bysort urn: egen q_old_fill = max(q_old)
bysort urn: egen q_new_fill = max(q_new)

tab q_old_fill q_new_fill
spearman q_old_fill q_new_fill
* Spearman rho = 0.77 - supports use of combined standardised measure

xtset urn year

* Year-on-year rank stability check
preserve

drop if year == 2020 | year == 2021
gen year_orig = year

* Recode years to fill Covid gap for lag construction
replace year = 2020 if year == 2022
replace year = 2021 if year == 2023
replace year = 2022 if year == 2024

gen q_year = .

levelsof year, local(years)
foreach y of local years {
    if `y' <= 2015 {
        xtile q_temp = ttaps if year == `y', nquantiles(5)
        replace q_year = q_temp if year == `y'
        drop q_temp
    }
    else {
        xtile q_temp = att8scr if year == `y', nquantiles(5)
        replace q_year = q_temp if year == `y'
        drop q_temp
    }
}

sort urn year
by urn: gen q_lag = q_year[_n-1]

sort urn year

levelsof year, local(years)
foreach y of local years {
    quietly count if !missing(q_year) & !missing(q_lag) & year == `y'
    if r(N) > 0 {
        quietly spearman q_year q_lag if year == `y'
        display "Year `y': rho = " r(rho) " (n = " r(N) ")"
    }
    else {
        display "Year `y': no observations"
    }
}

replace year = year_orig
drop year_orig q_year q_lag

restore

xtset urn year
capture drop mean_ttaps sd_ttaps z_ttaps mean_att8scr sd_att8scr z_att8scr

* Result: Year-on-year rank correlation ranges from 0.75-0.88
* 2016 transition within normal range - supports use of combined measure


********************************************************************************
**#                   3. Construction of Treatment Variables                  **
********************************************************************************

xtset urn year

* Binary treatment: library count increase/decrease within 1km
bysort urn (year): gen lib_first_1km = n_lib_1km if _n == 1
bysort urn: egen baseline_libs_1km  = max(lib_first_1km)
gen lib_change1 = n_lib_1km - baseline_libs_1km
drop lib_first_1km baseline_libs_1km

gen increase1 = (lib_change1 > 0)
gen decrease1 = (lib_change1 < 0)

* 2km radius (robustness)
bysort urn (year): gen lib_first_2km = n_lib_2km if _n == 1
bysort urn: egen baseline_libs_2km  = max(lib_first_2km)
gen lib_change2 = n_lib_2km - baseline_libs_2km
drop lib_first_2km baseline_libs_2km

gen increase2 = (lib_change2 > 0)
gen decrease2 = (lib_change2 < 0)

* Continuous distance treatment: change in log distance to nearest library
bysort urn (year): gen first_dist = dist_nearest_lib_km if _n == 1
bysort urn: egen baseline_dist    = max(first_dist)
gen change_dist = dist_nearest_lib_km - baseline_dist
drop first_dist

gen log_dist = log(dist_nearest_lib_km + 0.01)
bysort urn (year): gen log_dist_first = log_dist if _n == 1
bysort urn: egen log_baseline         = max(log_dist_first)
gen change_log_dist = log_dist - log_baseline
drop log_dist_first

* Binary distance increase/decrease
gen dist_increase = (change_dist > 0)
gen dist_decrease = (change_dist < 0)

* Strict distance treatment threshold (>0.5 log km decrease)
gen dist_treated_strict = (change_log_dist < -0.5)
gen dist_never          = (change_log_dist >= 0 | missing(change_log_dist))







********************************************************************************
**#                   4. Additional Controls and Flags                        **
********************************************************************************

* Parse school open and close dates from string
replace opendate  = "" if opendate  == "NA"
replace closedate = "" if closedate == "NA"

gen open_year  = .
gen close_year = .

replace open_year  = real(substr(opendate,  7, 4)) if regexm(opendate,  "^[0-9]{2}-[0-9]{2}-[0-9]{4}$")
replace close_year = real(substr(closedate, 7, 4)) if regexm(closedate, "^[0-9]{2}-[0-9]{2}-[0-9]{4}$")
replace open_year  = real(substr(opendate,  7, 4)) if regexm(opendate,  "^[0-9]{2}/[0-9]{2}/[0-9]{4}$")
replace close_year = real(substr(closedate, 7, 4)) if regexm(closedate, "^[0-9]{2}/[0-9]{2}/[0-9]{4}$")
replace open_year  = real(substr(opendate,  1, 4)) if regexm(opendate,  "^[0-9]{8}$")
replace close_year = real(substr(closedate, 1, 4)) if regexm(closedate, "^[0-9]{8}$")




* Drop observations outside school operating window
drop if open_year  != . & year < open_year
drop if close_year != . & year > close_year




* Flag schools that moved location (postcode change)
preserve
bysort urn (year): keep if _n == 1
keep urn postcode
rename postcode initial_postcode
tempfile postcodes
save `postcodes'
restore

merge m:1 urn using `postcodes', keep(master match) nogen

gen postcode_changed = (postcode != initial_postcode) & !missing(postcode)
bysort urn: egen school_moved = max(postcode_changed)
drop postcode_changed

* Nonlinear disadvantage variables
gen ptfsm6cla1a_sq = ptfsm6cla1a^2
xtile fsm_q = ptfsm6cla1a, nquantiles(5)

* Breakeven FSM share
gen break_even = (ptfsm6cla1a >= 0.23)

* Period indicator (pre/post score measure change)
gen post = (year > 2015)

* First treatment year for Callaway-Sant'Anna estimator
bysort urn (year): gen temp_cs = year if increase1 == 1
bysort urn: egen cs_treat_year = min(temp_cs)
replace cs_treat_year = 0 if missing(cs_treat_year)
drop temp_cs

* Baseline school score
bysort urn (year): gen baseline_score = overall_std if _n == 1
bysort urn: egen school_baseline = max(baseline_score)

xtset urn year

* Squared LA earnings for nonlinear controls
gen la_earnings_sq = la_earnings^2





********************************************************************************
**#               5. Main Regressions - Library Increase (1km)               **
********************************************************************************

* 5.1 Baseline - no disadvantage interaction
* Average effect is insignificant: masks heterogeneity by school disadvantage
reghdfe overall_std increase1 ///
    totpups ptfsm6cla1a ptealgrp2, ///
    absorb(urn year#la) vce(cluster urn)

* 5.2 Preferred specification: disadvantage interaction
* School and LA-year fixed effects; standard errors clustered at school level
reghdfe overall_std increase1##c.ptfsm6cla1a ///
    totpups ptfsm6cla1a ptealgrp2, ///
    absorb(urn year#la) vce(cluster urn)

* Calculate breakeven FSM share at which treatment effect changes sign
local treat_coef    = _b[1.increase1]
local interact_coef = _b[1.increase1#c.ptfsm6cla1a]
local breakeven     = `treat_coef' / abs(`interact_coef')
display "Breakeven FSM share: " `breakeven'
* Result: Increase coef 0.23***, interaction -0.57**, breakeven FSM ~0.40
* Library openings benefit lower-disadvantage schools; effect reverses above breakeven

* 5.3 By period (pre/post attainment measure change)
reghdfe overall_std increase1##c.ptfsm6cla1a ///
    totpups ptfsm6cla1a ptealgrp2 if year <= 2015, ///
    absorb(urn year#la) vce(cluster urn)

reghdfe overall_std increase1##c.ptfsm6cla1a ///
    totpups ptfsm6cla1a ptealgrp2 if year > 2015, ///
    absorb(urn year#la) vce(cluster urn)
* Results significant and comparable in each period

* 5.4 Test coefficient equality across periods
reghdfe overall_std increase1##c.ptfsm6cla1a ///
    increase1#i.post ///
    increase1#c.ptfsm6cla1a#i.post ///
    totpups ptfsm6cla1a ptealgrp2 post, ///
    absorb(urn year#la) vce(cluster urn)
* Cannot reject equality across periods

* 5.5 By disadvantage group (above/below breakeven FSM share)
reghdfe overall_std increase1##c.ptfsm6cla1a ///
    totpups ptfsm6cla1a ptealgrp2 if break_even == 0, ///
    absorb(urn year#la) vce(cluster urn)

reghdfe overall_std increase1##c.ptfsm6cla1a ///
    totpups ptfsm6cla1a ptealgrp2 if break_even == 1, ///
    absorb(urn year#la) vce(cluster urn)

* 5.6 School-specific trends robustness check
reghdfe overall_std increase1##c.ptfsm6cla1a ///
    totpups ptfsm6cla1a ptealgrp2, ///
    absorb(urn year#la c.year#urn) vce(cluster urn)
* Result strengthens to 0.252***, interaction -0.647***
* Confirms results not driven by pre-existing school-level trends
