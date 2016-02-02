*Dev stata 14
*Bryce Bartlett
*Robust estimates from time series data
*I've played around with this a lot, and there is no diff with OLS
*Need to figure out how to present by providing a test or two in the text
*alhtough, haven't really looked at interactions with IC

*@@@@@@@@@@@@@@
*Load Data (cleanded in R
*@@@@@@@@@@@@@@

*empty cache
clear
eststo clear

*@@@@@@
*data for cross-sectional asa table
*@@@@@@

import delimited ///
using "H:/projects/mort1/output/stata-series.csv", ///
numericcols(_all)

*egen icd10 cut(years), at 0 icodes

gen icd10=0
replace icd10 = 1 if(years>-1)

gen fyears = years + 6

*fe years

eststo: reg yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 i.fyears female complex home ltcare oplace black

eststo: reg yrrdc a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 i.fyears female complex home ltcare oplace black

esttab using "H:/projects/mort1/output/fixed-years.rtf", ///
 scalars(rank ll) replace wide b(3) se(3) t
 
eststo clear
 
*fe demean style across cells

xtset cell

eststo: xtreg yrrac years, fe

eststo: xtreg yrrac years icd10 c.years#c.icd10, fe

eststo: xtreg yrrdc years, fe

eststo: xtreg yrrdc years icd10 c.years#c.icd10, fe

esttab using "H:/projects/mort1/output/fixed-cells.rtf", ///
 scalars(rank ll) replace wide b(3) se(3) t
 
eststo clear


*@@@@@@
*cross-classified model 1
*@@@@@@

*as is the age is fixed...
*eststo clear
*eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
* years female complex home ltcare oplace black ///
* || _all: R.a || _all: R.years, var reml

*eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
* years female complex home ltcare oplace black ///
* || cell:, var reml

*eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
* icd10 c.icd10#(a45 a50 a55 a60 a65 a70 a75 a80 a85) ///
* years icd10#c.years female c.icd10#female black c.icd10#black ///
* complex home ltcare oplace || cell:, var reml

*two varance components:
*year is very small (0.005 with cmall covariance (similar scale))

eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 years female complex home ltcare oplace black ///
 || cell: years, var cov(un) reml

eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 icd10 c.icd10#(c.a45 c.a50 c.a55 c.a60 c.a65 c.a70 c.a75 c.a80 c.a85) ///
 years c.icd10#c.years female c.icd10#c.female black c.icd10#c.black ///
 complex home ltcare oplace || cell: years, var cov(un) reml

eststo: xtmixed yrrdc a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 years female complex home ltcare oplace black ///
 || cell: years, var cov(un) reml

eststo: xtmixed yrrdc a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 icd10 c.icd10#(c.a45 c.a50 c.a55 c.a60 c.a65 c.a70 c.a75 c.a80 c.a85) ///
 years c.icd10#c.years female c.icd10#c.female black c.icd10#c.black ///
 complex home ltcare oplace || cell: years, var cov(un) reml

esttab using "H:/projects/mort1/output/hierarchical.rtf", ///
 scalars(rank ll) replace wide b(3) se(3) t
 
eststo clear
