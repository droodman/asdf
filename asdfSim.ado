*! adsf 0.1.1 23 May 2020
*! Copyright (C) 2020 David Roodman

// asdf: Analytical Stochastic Dynamic Framework
// Routines for fitting and simulating univariate stochastic models with analytical statements for their distributions--power of CIR/Feller, geometric Brownian motion
// Also does "superexponential" NLS.

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

* generate simulated quantiles, sample paths, and final distribution from current estimate. Leave behind in Mata asdfSim* variables.
cap program drop asdfSim
program define asdfSim
	syntax, t(real) tsamp(integer) tres(integer) m(integer) y0(real) t0(real) [nq(integer 0) ns(integer 0) nov model(string) obsquant]
	mata tmp = asdfSimPaths(`t', `tsamp', `tres', `m', `y0', `t0', `nq', `ns', ("`v'"!=""? 0 : st_matrix("e(V)")), asdfSim`e(model)'(), "`obsquant'"!="")
  mata asdfSimt = *tmp[1]; asdfSimY = *tmp[2]; asdfSimYf = *tmp[3]; asdfSimYq = *tmp[4]; asdfSimYs = *tmp[5]; asdfSimYqObs = *tmp[6]
end
