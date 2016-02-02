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


gen fyears = years + 6

*fe years

eststo: reg yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 i.fyears female complex home ltcare oplace black

*fe demean style with a45

xtset a

xtreg yrrac years, fe

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
 || _all: R.years, var reml
