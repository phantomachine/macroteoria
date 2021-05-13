clear
set more off
set mem 400m
set logtype text
capture log close

cd C:\temp

#delimit ;
clear;

use Data;

gen tight =0;
replace uirate = uirate / 100;
replace tight = vacancies/uirate;

sort year quarter;
drop if year < 1951;
drop if year > 2004;

gen wage_gen=labor_share*output_pp;
gen log_wage_gen=log(wage_gen);
gen log_output_pp = log(output_pp);
gen ui_q= uirate;
gen v_q= vacancies;
gen tight_q= tight;

gen time=4*(year-51)+quarter;
tsset time;

******************************;
**** Wage Regression   *******;
******************************;
hprescott log_wage_gen, stub(H) smooth(1600);
hprescott log_output_pp, stub(H) smooth(1600);
reg H_log_wage_gen_1 H_log_output_pp_1;
******************************;

gen log_output_per_p = 0;
gen filter_output_p = 0;
replace log_output_per_p = log(output_pp);
hprescott log_output_per_p, stub(H) smooth(1600);
replace filter_output_p = H_log_output_per_p_1;

gen log_ui_q = log(ui_q);
gen filter_ui_q = 0;
hprescott log_ui_q, stub(H) smooth(1600);
replace filter_ui_q = H_log_ui_q_1;

gen log_v_q = log(v_q);
gen filter_v_q = 0;
hprescott log_v_q, stub(H) smooth(1600);
replace filter_v_q = H_log_v_q_1;

gen log_tight_q=0;
replace log_tight_q = log(tight_q);
gen filter_tight_q = 0;
hprescott log_tight_q, stub(H) smooth(1600);
replace filter_tight_q = H_log_tight_q_1;

******************************;
**** Standard Deviations *****;
******************************;
summarize filter_v_q;
summarize filter_ui_q;
summarize  filter_tight_q;
summarize  filter_output_p;

***************************;
*** Correlation Matrix ****;
***************************;
correlate filter_ui_q filter_v_q filter_tight_q filter_output_p;

*************************;
*** Autocorrelation  ****;
*************************;
gen L_find_q = 0;
gen L_u_q = 0;
gen L_v_q = 0;
gen L_tight_q = 0;
gen L_pp_q = 0;

replace L_u_q = filter_ui_q[_n+1];
replace L_v_q = filter_v_q[_n+1];
replace L_pp_q= filter_output_p[_n+1];
replace L_tight_q = filter_tight_q[_n+1];

drop if year==2004 & quarter==4;

correlate L_u_q filter_ui_q;
correlate L_v_q filter_v_q;
correlate L_tight_q filter_tight_q;
correlate L_pp_q filter_output_p;

