*! adsf 0.1 23 May 2020
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

cap program drop asdf_p
program define asdf_p
  syntax anything, [cdf *]

  if "`cdf'"=="" {
  	ml_p `anything', `options'
    exit
  }

  tempvar x0 tdelta prevobs hold
 	_estimates hold `hold', restore copy  // surface the sample marker without duplicating it(?)
  qui if "`e(modeltype)'" == "dynamic" {
		gen long `prevobs' = . in 1
		replace `prevobs' = cond(`e(tvar)'[_n-1]<. & `e(depvar)'[_n-1] <., _n-1, `prevobs'[_n-1]) if _n>1
		gen double `x0' = `e(depvar)'[`prevobs'] if `hold'
		gen double `tdelta' = `e(tvar)' - `e(tvar)'[`prevobs'] if `hold'
	}
  else if "`model'" != "bernounls" gen double `tdelta' = `e(tvar)' - `e(t0)' if `hold'

  mata S = asdfEst`e(model)'()
  mata S.setData("`e(modeltype)'"=="dynamic"? st_data(., "`x0'", "`hold'") : st_numscalar("`e(x0)'"), st_data(., "`e(depvar)'" , "`hold'"), st_data(., "`tdelta'" , "`hold'"))
  mata t=S.getParamFixedNames(); for(i=length(t);i;i--) S.setFixedParam(t[i], st_numscalar("e("+t[i]+")"))
  mata S.processParams(st_matrix("e(b)"))

  _score_spec `anything'
  qui gen double `anything' = .
  mata st_view(asdf_cdf=., ., "`s(varlist)'", "`hold'"); asdf_cdf[,] = S.CDF()
end