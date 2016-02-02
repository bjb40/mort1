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

*fe demean style across cells

xtset cell

xtreg yrrac years, fe

xtreg yrrac years icd10 c.years#c.icd10, fe

*@@@@@@
*cross-classified model 1
*@@@@@@

*as is the age is fixed...
*eststo clear
*eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
* years female complex home ltcare oplace black ///
* || _all: R.a || _all: R.years, var reml

eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 years female complex home ltcare oplace black ///
 || cell:, var reml

eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 icd10 c.icd10#(a45 a50 a55 a60 a65 a70 a75 a80 a85) ///
 years icd10#c.years female c.icd10#female black c.icd10#black ///
 complex home ltcare oplace || cell:, var reml

*two years

eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 years female complex home ltcare oplace black ///
 || cell: years, var reml

eststo: xtmixed yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 icd10 c.icd10#(a45 a50 a55 a60 a65 a70 a75 a80 a85) ///
 years icd10#c.years female c.icd10#female black c.icd10#black ///
 complex home ltcare oplace || cell: years, var reml
