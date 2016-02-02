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
using  "H:/projects/mort1/output/stata-series.csv", ///
numericcols(_all)

*@@@@@
*Analysis
*newey doesn't work b/c this is panel data
*this is just an iid problem
*@@@@@


xtset cell years

gen icd10 = 0
replace icd10=1 if years >= 0

*@@@@@
*Testing serial correlation
*have to refix the data to include only demographic variables
*@@@@@@

 *testing for unit root
 xtunitroot fisher yrrac, dfuller lags(1)
 
 *no evidence of unit root for yrrac; there is for yrrdc
 xtunitroot fisher yrrdc, dfuller lags(1)
 
*@@@@@
*Model 1,2 yrrac
*@@@@@

 
 *this is the fgls n woolrdge p. 421 / inefficient, and employes a
eststo: xtgls yrrac  a45 a50 a55 a60 a65 a70 a75 a80 a85 years ///
female complex home ltcare oplace black, panels(hetero) corr(psar1)

eststo: xtgls yrrac a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 icd10 c.icd10#(a45 a50 a55 a60 a65 a70 a75 a80 a85) ///
 years icd10#c.years female c.icd10#female black c.icd10#black ///
	, panels(hetero) corr(psar1)
 
 *delete structural zeros
 drop if missing(yrrdc)
 
 
*@@@@@
*Model 3,4 yrrdc
*@@@@@

eststo: xtgls yrrdc  a45 a50 a55 a60 a65 a70 a75 a80 a85 years ///
female complex home ltcare oplace black, panels(hetero) corr(psar1) force

eststo: xtgls yrrdc a45 a50 a55 a60 a65 a70 a75 a80 a85 ///
 icd10 c.icd10#(a45 a50 a55 a60 a65 a70 a75 a80 a85) ///
 years icd10#c.years female c.icd10#female black c.icd10#black ///
	, panels(hetero) corr(psar1) force
	
*@@@@@
*Output table
*@@@@@
	
esttab using "H:/projects/mort1/output/flgs-tables.rtf", ///
 scalars(rank ll) replace wide b(3) se(3) t

 
