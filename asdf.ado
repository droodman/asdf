*! adsf 0.1.1 23 May 2020
*! Copyright (C) 2020 David Roodman

// asdf: Analytical Stochastic Dynamic Framework
// Routines for fitting and simulating univariate stochastic models with analytical statements for their distributions--Bernoulli, geometric Brownian motion
// Also does Bernoulli NLS.

* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.

*! Version history at bottom

cap program drop asdf
program define asdf, eclass
	if replay() {
		if "`e(cmd)'" != "asdf" error 301
		if _by() error 190
		ml di
		exit 0
	}

	syntax varlist(min=1 max=2 numeric) [pw fw aw iw/] [if] [in], model(string) [MODELType(string) CONSTraints(passthru) TECHnique(passthru) from(namelist min=1 max=1) lf(string) CLuster(passthru) Robust vce(passthru) *]
  if "`modeltype'"=="" local modeltype dynamic
  if !inlist("`modeltype'", "static", "dynamic") error 198

	tokenize `varlist'
	local x `1'
	local tvar `2'

	tempvar touse
	qui mark `touse' `if' `in'
	qui if "`tvar'" == "" {  // if time var not provided, take it from tsset
		local tvar: char _dta[tis]
		if "`tvar'"=="" error 111
  }
	qui markout `touse' `tvar' `x'

	if "`weight'" != "" {
		cap confirm var `exp'
		if _rc {
			tempname wvar
			qui gen double `wvar' = `exp' if `touse'
		}
		else local wvar `exp'
		local wgt  [`weight'=`wvar']
		local wtype `weight'
		local wexp `"=`exp'"'
    qui markout `touse' `wvar'
	}

  tempname ti tf depvari depvarf
  sum `tvar' if `touse', meanonly
  scalar `ti' = r(min)
  scalar `tf' = r(max)
  sum `x' if `touse', meanonly
  scalar `depvari' = r(min)
  scalar `depvarf' = r(max)

  tempvar x0 tdelta prevobs
  qui if "`modeltype'" == "dynamic" {
		gen long `prevobs' = . in 1
		replace `prevobs' = cond(`tvar'[_n-1]<. & `x'[_n-1] <., _n-1, `prevobs'[_n-1]) if _n>1
		gen double `x0' = `x'[`prevobs'] if `touse'
		gen double `tdelta' = `tvar' - `tvar'[`prevobs'] if `touse'
    markout `touse' `x0'
	}
  else {
  	tempvar n t0
    gen long `n' = _n
    sum `n' if `touse', meanonly
    scalar `x0' = `x'[r(min)]
    scalar `t0' = `tvar'[r(min)]
    if "`model'" != "bernounls" gen double `tdelta' = `tvar' - `t0' if `touse' & _n > r(min)
  }
  markout `touse' `tdelta' `wvar'

	local _options `options'
	local 0, `model'
	syntax, [bernoudiff bernoudiff2 gbm bernounls stickyfeller]
	local options `_options'

	mata S = asdfEst`model'()
	if `"`lf'"'=="" mata st_local("lf", strofreal(S.getlf()))

  if "`model'" == "bernounls" {
    mata S.setData("`modeltype'"=="dynamic"? st_data(., "`x0'", "`touse'") : st_numscalar("`x0'"), st_data(., "`x'" , "`touse'"), st_data(., "`tdelta'" , "`touse'") `=cond("`wtype'" != "", `", st_data(., "`wvar'" , "`touse'")"', "")')
    tempvar b V
    if "`from'" != "" mata S.Estimate(st_matrix("`from'"))
                 else mata S.Estimate()
    mata st_matrix("`b'", S.b); st_matrix("`V'", S.V)
    if `b'[1,1] != . {
    	mata colnames = J(4,1,"/"),S.getParamEstNames()'; st_matrixcolstripe("`b'", colnames); st_matrixcolstripe("`V'", colnames); st_matrixrowstripe("`V'", colnames)
      qui count if `touse'
      ereturn post `b' `V', obs(`r(N)') esample(`touse')
      ereturn local converged = 1
      eret di
    }
    else error 430
  }
  else {
    mata S.setData("`modeltype'"=="dynamic"? st_data(., "`x0'", "`touse'") : st_numscalar("`x0'"), st_data(., "`x'" , "`touse'"), st_data(., "`tdelta'" , "`touse'"))
    mata st_local("params", invtokens(S.getParamEstNames()))
    foreach param in `params' {
      local syntaxopts `syntaxopts' `=upper("`param'")'VARs(varlist ts fv)
    }
    local 0, `options'
    syntax, [`syntaxopts' *]

    foreach param in `params' {
      markout `touse' ``param'vars'
      local cmdline `cmdline' `=cond("``param'vars'"=="", "/`param'", "(`param': ``param'vars')")'
    }
    
    mata st_local("fixedparams", invtokens(S.getParamFixedNames()))
    foreach param in `fixedparams' {
      local 0, `options'
      syntax, `param'(real) [*]
      mata S.setFixedParam("`param'", ``param'')
    }

    ml model lf`lf' asdflf2() `cmdline' if `touse' `wgt', `constraints' `technique' nopreserve `robust' `cluster' `vce'
    mata moptimize_init_userinfo($ML_M, 1, S)
    mata moptimize_init_nmsimplexdeltas($ML_M, 1e-3)  // ? in case tech(nm)
    if "`from'" != "" ml init `from', copy
    mlopts mlopts, `options'
    ml max, `mlopts' search(off)
  }
  
	foreach retval in `fixedparams' ti tf depvari depvarf {
		ereturn scalar `retval' = ``retval''
	}
	foreach param in `params' {
		ereturn local `param'vars ``param'vars'
	}
	ereturn local model `model'
	ereturn local modeltype `modeltype'
  ereturn local depvar `x'
  ereturn local tvar `tvar'
  if "`modeltype'" == "static" {
    ereturn scalar x0 = `x0'
    ereturn scalar t0 = `t0'
  }
  ereturn local wexp `wexp'
  ereturn local predict asdf_p
	ereturn local cmd asdf
end
