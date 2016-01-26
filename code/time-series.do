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
using "H:/Academic Projects/Mortality I/output/stata-series.csv", ///
numericcols(_all)

*@@@@@
*Analysis
*newey doesn't work b/c this is panel data
*this is just an iid problem
*@@@@@


xtset dem years

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
eststo: xtgls yrrac v6-v14 years female complex home ltcare oplace black, panels(hetero) corr(psar1)

eststo: xtgls yrrac v6-v14 years female complex home ltcare oplace black ///
	icd10 c.icd10#(v6-v14) c.icd10#c.years ///
	, panels(hetero) corr(psar1)
 
 *delete structural zeros
 drop if missing(yrrdc)
 
 
*@@@@@
*Model 3,4 yrrdc
*@@@@@

eststo: xtgls yrrdc v6-v14 years female complex home ltcare oplace black, panels(correlated) force

eststo: xtgls yrrdc v6-v14 years female complex home ltcare oplace black ///
	icd10 c.icd10#(v6-v14) c.icd10#c.years ///
	, panels(hetero) force
	
*@@@@@
*Output table
*@@@@@
	
esttab using "H:/Academic Projects/Mortality I/output/flgs-tables.rtf", ///
 scalars(rank ll) replace wide b(3) se(3) t

 
*@@@@@
*Reload for general (this is more like an HAPC model...
*BUT you don't have to go exaclty there...
*@@@@@
  
clear
eststo clear

*@@@@@@
*data for cross-sectional asa table
*@@@@@@

import delimited ///
using "H:/Academic Projects/Mortality I/output/stata-series.csv", ///
numericcols(_all)

*@@@@@
*Analysis
*Nested Hierarchical and cross-classified
*@@@@@


*xtset dem years

gen icd10 = 0
replace icd10=1 if years >= 0

xtset, clear

*@@@@@
*Models 1 and 2
*@@@@@
eststo clear
eststo: xtmixed yrrac v6-v14 years female complex home ltcare oplace black ///
 || _all: R.dem || _all: R.years, var reml
 
eststo: xtmixed yrrac v6-v14 c.icd10#(v6-v14) years female complex home ltcare oplace black ///
  icd10 c.icd10#c.years ///
 || _all: R.dem || _all: R.years, var reml
	

*@@@@@
*Models 3 and 4
*@@@@@

eststo: xtmixed yrrdc v6-v14 years female complex home ltcare oplace black ///
  || _all: R.dem || _all: R.years, var reml
 
eststo: xtmixed yrrdc v6-v14 c.icd10#(v6-v14) years female complex home ltcare oplace black ///
  icd10 c.icd10#c.years ///
   || _all: R.dem || _all: R.years, var reml


*@@@@@
*output results
*@@@@@

esttab using "H:/Academic Projects/Mortality I/output/ccrem-tables.rtf", ///
 scalars(rank ll) replace wide b(3) se(3) t
	
*@@@@@
*Get margins for age-specific
*@@@@@



