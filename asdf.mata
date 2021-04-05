*! asdf 0.1 23 May 2020
*! Copyright (C) 2020 David Roodman

// asdf: Analytical Stochastic Dynamic Framework
// Routines for fitting and simulating univariate stochastic models with analytical statements for their distributions--Bernoulli, geometric Brownian motion
// Also does Bernoulli NLS.

// Requires Stata 16 or later

* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.

//
// Class asdfSim is for simulation, asdfEst for estimation. Both are extended to child classes specific to each model.
//

// Compile-time parameters for sticky Feller routines
global BesselKPower_M 22 // max term index in power series expansion of BesselK()
global BesselKAsym_M 7   // max term index in asymptotic approximation of BesselK()
global UPower_M 60       // max term index in power series expansion of U()
global UAbadSesma_N 20   // number of terms in Abad-Sesma large-parameter approximation of U()

scalar log_epsilon = ln(1e-15)                        // Logic governing evaluation points on parabola of integration, extracted from Garrappa's mlf routine for Matlab
mata st_numscalar("log_eps", ln(epsilon(1)))
scalar _mu = log_epsilon - log_eps                    // mu as in Trefethen and Weideman, scaling parabola, not stickiness parameter mu
scalar w = sqrt(-log_eps / _mu)                       // half-width of integration range, needed to assure given precision
scalar N = ceil(-w * log_epsilon / (2*_pi)) + 1       // number of integration points (on positive arm of parabola)
scalar h = w / (N-.5)                                 // width of bars in Riemann integral; do an even number to avoid double-counting evaluation at vertex
global StickyFeller_h  = substr(strofreal(  h, "%21x"), 1 + (  h >= 0), .)  // strip off any leading "+" to prevent syntax errors
global StickyFeller_mu = substr(strofreal(_mu, "%21x"), 1 + (_mu >= 0), .)
global StickyFeller_N = N

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

struct smatrix {
	real matrix M
}


// general SDE simulator class
class asdfSim {
	real scalar reflect, dirty
	real colvector Yf, t  // // final distribution; time values for simulated paths and quantiles
	real matrix Y, Yq  // sample paths, and quantiles of simulated paths
	real rowvector Yi  // temp var: current cross-section of simulations
	real colvector b
	real matrix V
	real scalar Y0, t0, T, Tsamp, Tres, M

	private void sim()
	void new(), setSimParams(), setModelParams()
	virtual void simStart(), simStep(), simStop()
	real matrix Y(), Yq()
	real colvector Yf(), t(), YqObs()
}

void asdfSim::new() dirty = 1

void asdfSim::setModelParams(real colvector _b, real matrix _V, | real scalar _reflect) {
	b = _b; V = _V; reflect = _reflect
	dirty = 1
}

void asdfSim::setSimParams(real scalar _Y0,  // starting value
								           real scalar _t0,  // starting time--only affects labelling of time coordinates
							             real scalar _T,  // time span of simulation
							             real scalar _Tsamp,  // number of evenly spaced time points to retain and return; 1 = just endpoints 
							             real scalar _Tres,  // time steps between retained points
							             real scalar _M) {  // number of paths to simulate
	Y0 = _Y0; t0 = _t0; T = _T; Tsamp = _Tsamp; Tres = _Tres; M = _M
	dirty = 1
}

// run a simulation
void asdfSim::sim() {
	real rowvector negative; real scalar i, j, s; real matrix X; real rowvector L
  pragma unset X; pragma unset L

	Yi = J(1, M, Y0)
	i = Tres * Tsamp

	if (V==. | V==0)
		simStart(b, T/i)
	else {
		symeigensystem(V, X, L)
		simStart(b :+ (X :* sqrt(edittozero(L,1))) * rnormal(rows(b), M, 0, 1), T/i)  // X :* sqrt(L) is sqrt(V) computed in way robust to V being singular, when estimation is constrained
	}

	Y  = J(Tsamp+(Tsamp > 1), 1, Yi)  // start each path at starting point, unless only final values requested

	j = 1 + (Tsamp > 1)
	s = Tres - 1

	for (; i; i--) {
		simStep(reflect)
		negative = selectindex(Yi :< 0)
		if (length(negative)) Yi[negative] = reflect? -Yi[negative] : J(1, length(negative), 0)  // reflect (sort of) or absorb boundary crossers; reflection method is noticeably problematic as nu->-1 and theoretical density at zero goes to infinity
		if (--s <= 0) {
			Y[j++,] = Yi
			s = s + Tres
		}
	}

	simStop()
	t = rangen(t0, t0+T, Tsamp+1)  // time counter
	dirty = 0
}

// return full set of simulated paths
real matrix asdfSim::Y() {
	if (dirty) sim()
	return (Y)
}

// evenly spaced quantiles of simulated paths. E.g., Nq = 20 --> 5th, 10th, etc.
real matrix asdfSim::Yq(real scalar Nq) {
	real scalar i
	if (dirty) sim()
	if (Nq<. & Nq>0) {
		Yq = Y'  // quantiles of distribution at each time step, treating missing as +infty
		for (i=cols(Yq);i;i--)
			Yq[,i] = sort(Yq[,i],1)
		Yq = Yq[round((1..Nq-1)/Nq * M),]'
		return (Yq)
	}
}

// quantiles of specific time/value points, with linear interpolation between sampling points
real colvector asdfSim::YqObs(real colvector tObs, real colvector YObs) {
	real colvector _t, floor, frac
	if (dirty) sim()
	if (rows(tObs)) {
    _t = (tObs :- t0) * (Tsamp / T) :+ 1
    floor = floor(_t)
    frac = _t - floor
    return ((rowsum(Y[floor,] :< YObs) :* (1 :- frac) + rowsum(Y[floor:+1,] :< YObs) :* frac) / M)  // will crash if an element of tObs == t0 + T (final time point)
	}
}

// final values of all paths
real colvector asdfSim::Yf() {
	if (dirty) sim()
	return (Y[rows(Y),]')
}

// simulation sampling time points
real colvector asdfSim::t() {
	if (dirty) sim()
	return (t)
}


// Child class for simulating Bernoulli SDE by simulating Feller diffusion dX = (bX+c) + sqrt(2aX)dW then taking Y = X^gamma
class asdfSimbernoudiff extends asdfSim {
	real rowvector _a, _b, _c, gamma  // temp vars
	real scalar colsYicolsa
	virtual void simStart(), simStep(), simStop()
}

// model-specific prep steps for simulation
void asdfSimbernoudiff::simStart(real matrix params, real scalar dt) {
	_c =      dt * (1 :+ params[3,]) :* (_a = exp(params[1,]))  // c = (1+nu)*a
	_a = sqrt(dt * (_a + _a))  // _a
	_b = 1 :+ dt * params[2,]  // _b = 1 + b
	gamma = params[4,]
	Yi = Yi :^ (1 :/ gamma)
	colsYicolsa = cols(Yi)/cols(_a)
}

// model-specific simulation step
void asdfSimbernoudiff::simStep(real scalar reflect)
	Yi = _b :* Yi :+ (reflect? _c : _c :* (Yi:!=0)) :+ sqrt(Yi) :* rnormal(1, colsYicolsa, 0, _a)

// model-specific finalization
void asdfSimbernoudiff::simStop()
	Y = Y :^ gamma


// Class for simulating Bernoulli SDE via s, B, delta, lnsig parameterization

class asdfSimbernoudiff2 extends asdfSimbernoudiff {
	virtual void simStart()
}

// model-specific prep steps for simulation
void asdfSimbernoudiff2::simStart(real matrix params, real scalar dt) {
	real rowvector t, s, B, delta, lnsig, nu

	s = params[1,]; B = params[2,]; delta = params[3,]; lnsig = params[4,]
	t = lnsig + ln(abs(B))
	_a = exp(t + t  :- 1.62e42fefa39efX-001 /*ln(2)*/)
	_b = -B :* delta
	nu = (1 :- (s+s) :/ exp(lnsig+lnsig)) :/ B
	gamma = -1 :/ B

	_c =      dt * (1 :+ nu) :* _a  // c = (1+nu)*a
	_a = sqrt(dt * (_a + _a))  // _a
	_b = 1 :+ dt * _b  // _b = 1 + b
	Yi = Yi :^ (-B)
	colsYicolsa = cols(Yi)/cols(_a)
}


// Class for simulating geometric Brownian motion
class asdfSimgbm extends asdfSim {
	real rowvector _a, _b  // temp vars
	virtual void simStart(), simStep(), simStop()
}

// model-specific prep steps for simulation
void asdfSimgbm::simStart(real matrix params, real scalar dt) {
	_a = exp(params[1,]); _a  = sqrt(dt * (_a+_a))
	_b = 1 :+ dt * params[2,]
}

// model-specific simulation step
void asdfSimgbm::simStep(| real scalar reflect) {  // argument not used
	pragma unset reflect
	Yi = (_b :+ rnormal(1, cols(Yi)/cols(_a), 0, _a)) :* Yi
}

void asdfSimgbm::simStop() {}

// do the work of post-estimation command asdfSim.ado
pointer (real matrix) colvector asdfSimPaths(real scalar T, real scalar Tsamp, real scalar Tres, real scalar M, real scalar Y0, real scalar T0, real scalar Nq, real scalar Ns, real matrix V, class asdfSim scalar S, real scalar obsquant) {
  real matrix Ys, YqObs
  real colvector _tObs, _YObs, tObs, YObs, samp
  pragma unset _tObs; pragma unset _YObs; pragma unset tObs; pragma unset YObs

  S.setModelParams(st_matrix("e(b)")', V, st_numscalar("e(reflect)"))
  S.setSimParams(Y0, T0, T, Tsamp, Tres, M)
  if (Ns > 0)
    Ys = S.Y()[|.,.\.,Ns|]
  if (obsquant) {
    st_view(_tObs, ., st_global("e(tvar)"))
    st_view(_YObs, ., st_global("e(depvar)"))
    samp = _tObs :< . :& _YObs :< .
    st_select(tObs, _tObs, samp)
    st_select(YObs, _YObs, samp)
    YqObs = S.YqObs(tObs, YObs)
  }
  return (&(S.t()) \ &(S.Y()) \ &(S.Yf()) \ &(S.Yq(Nq)) \ &Ys \ &YqObs)  // time, full path set, final path values, quantiles of full distribution, quantiles of provided observations
}


//
// Estimator classes
//

class asdfEst {  // SDE model with analytically expressible distribution
	real scalar reflect, Nzeros
	real colvector lnY0, lnY, tDelta, tDelta2, zeroind, zerofill
	struct smatrix colvector params

	virtual struct smatrix colvector getParams()
	real scalar getlf()
	void setData()
	virtual void lnPDF(), setFixedParam()
	virtual string rowvector getParamEstNames(), getParamFixedNames()
}

real scalar asdfEst::getlf() return(2)  // default is that estimator is type lf2

// _Y0, _Y, _tDelta = inital & final values, and time gaps
void asdfEst::setData(real colvector _tDelta, real colvector _Y0, real colvector _Y) {
	lnY0 = ln(_Y0)
	lnY  = ln(_Y )
	zerofill = J(Nzeros = length(zeroind = selectindex(!_Y)), 1, 0)  // not used in GBM
	tDelta = _tDelta; tDelta2 = tDelta :* tDelta
}

struct smatrix colvector asdfEst::getParams() return(params)


// Bernoulli NLS estimator class

class asdfEstbernounls extends asdfEst {
	real colvector Ydot, wtYdot
  pointer(real colvector) scalar pwt, pY
  real rowvector b
  real matrix V

	void Estimate(), d2(), setData()
  virtual string rowvector getParamEstNames()
}

string rowvector asdfEstbernounls::getParamEstNames()
	return (("s", "B", "delta", "lnsig"))

// _Y0, _Y, _tDelta = inital & final values, and time gaps; wt = optional weights
void asdfEstbernounls::setData(real colvector _tDelta, real colvector _Y0, real colvector _Y, | real colvector wt) {
	Ydot = (_Y :/ _Y0) :^ (1:/_tDelta) :- 1
	pY = &_Y0
  lnY = ln(_Y0)
	pwt = rows(wt)? &(_tDelta :* wt) : &_tDelta  // weighting by time gap as in Kremer (1993) along with any user-supplied weights
	wtYdot = *pwt :* Ydot
}


// essentially concentrated log likelihood for just B, in Bernoulli NLS; not really log-likelihood but corresponding sum of squared errors
void BernoulliNLSc(real scalar todo, real scalar B, class asdfEstbernounls S, real scalar sig2, real matrix g, real matrix H)
  S.d2(todo, B, sig2, g, H)

void asdfEstbernounls::d2(real scalar todo, real scalar B, real scalar sig2, real matrix g, real matrix H) {
	real colvector YB, e, ewt, Xdot, t2, Ht2, wtXdot; real scalar beta, ewtXdot
	YB = *pY :^ B
	beta = (H = invsym(quadcross(YB, 1, *pwt, YB, 1))) * quadcross(YB, 1, wtYdot, 0)
	e = Ydot :- YB * beta[1] :- beta[2]
	ewt = e :* *pwt
	sig2 = ewt ' e
	if (todo) {
		beta = beta[1]
		Xdot = YB :* lnY
    ewtXdot = quadcross(ewt, Xdot)
		g = -2 * beta * ewtXdot
		if (todo > 1) {
			wtXdot = *pwt :* Xdot
			t2 = quadcross(YB, 1, wtXdot, 0)
			Ht2 = H * t2
			H =  beta * (beta * (quadcross(Xdot, wtXdot) - quadcross(t2, Ht2)) - quadcross(ewt, Xdot, lnY)) - H[1,1]*ewtXdot*ewtXdot - g*Ht2[1]
			H = H + H
		}
	}
}

// weighted NLS estimatation of the Bernoulli diff eq
// if initial search doesn't converge, stores . in B and computes other results as if B=0
void asdfEstbernounls::Estimate(| real scalar initB) {
	transmorphic scalar S; real colvector YB, e, ds, dB, slnY, YB2, YB2sinvsig2; real scalar B, rc, sig2, sumwt

	S = optimize_init()
	optimize_init_evaluator(S, &BernoulliNLSc())
	optimize_init_params(S, editmissing(initB, 1))  // default starting value = 1
  optimize_init_argument(S, 1, this)
	optimize_init_which(S, "min")
	optimize_init_evaluatortype(S, "d2")
  if (rc = _optimize(S)) {
    b = .
    return
  }
	B = optimize_result_params(S)  // concentrated ML estimate of B
	YB = *pY :^ B
	b = lusolve(quadcross(YB, 1, *pwt, YB, 1), quadcross(YB, 1, *pwt, Ydot, 0))
	e = (Ydot - b[1] * YB) :- b[2]
  sumwt = rows(*pwt)>1 ? quadcolsum(*pwt) : rows(*pY)
  sig2 = quadcrossdev(e, mean(e, *pwt), *pwt, e, 0) / sumwt
	b = (b[1] , (rc? . : B) , b[2] , ln(sig2)*.5)

  V = J(4,4,0)
  dB = (ds = e :* YB) :* (slnY = b[1] * lnY)
  YB2 = YB :* YB; YB2sinvsig2 = YB2 * b[1]
  V[1,1] = quadcross(*pwt, YB2)
  V[1,2] = quadcross(*pwt, lnY:*(YB2sinvsig2 - ds))
  V[1,3] = quadcross(*pwt, YB)
  V[1,4] = quadcross(*pwt, ds) * 2
  V[2,2] = quadcross(*pwt, lnY:*(YB2sinvsig2 :* slnY - dB))
  V[2,3] = quadcross(*pwt, slnY:*YB)
  V[2,4] = quadcross(*pwt, dB) * 2
  V[3,3] = sumwt
  V[3,4] = quadcross(*pwt, e) * 2
  V[4,4] = quadcross(*pwt, e, e) * 2
  V = invsym(V) * (sig2 * sumwt / rows(*pY))
}


// class for estimation of Bernoulli diffusion with lna, b, nu, gamma parameterization

class asdfEstbernoudiff extends asdfEst {
	pointer (real colvector) scalar pb, plna, pnu, pgamma, psigma2
	real colvector bt, lnlambda, lnx, lnC, lnX0, lnX, dlnlambdadgamma, dlnxdgamma, bdivexpm1bt

	virtual void lnPDF(), setFixedParam(), getScalarParams(), getSmatrixParams(), transformGradient(), processParams()
	virtual string rowvector getParamEstNames(), getParamFixedNames()
	void new()
  real colvector CDF(), invCDF()
}

string rowvector asdfEstbernoudiff::getParamEstNames()
	return (("lna", "b", "nu", "gamma"))

string rowvector asdfEstbernoudiff::getParamFixedNames()
	return ("reflect")

void asdfEstbernoudiff::setFixedParam(string scalar name, numeric value) {
	if (name == "reflect") reflect = value
}

void asdfEstbernoudiff::new()
	params = smatrix(4,1)  // pre-allocated parameter vessel
	
void asdfEstbernoudiff::getScalarParams(real vector params) {
	plna = &(params[1]); pb = &(params[2]); pnu = &(params[3]); pgamma = &(params[4])
}

void asdfEstbernoudiff::getSmatrixParams(struct smatrix vector params) {
	plna = &(params[1].M); pb = &(params[2].M); pnu = &(params[3].M); pgamma = &(params[4].M)
}

void asdfEstbernoudiff::processParams(transmorphic vector params) {
	real colvector ind

	if (eltype(params)=="real")
		getScalarParams(params)
	else
		getSmatrixParams(params)

	bt = *pb :* tDelta
	if (*pb==0) {
		lnC = -ln(tDelta)
    bdivexpm1bt = 0
	} else {
		lnC = ln(bdivexpm1bt = editmissing(*pb :/ expm1(bt), 0))
		if (length(ind = selectindex(!*pb)))
			lnC[ind] = -ln(tDelta[ind])
		if (length(ind = selectindex(bt:>1.628b76e3a78f2X+009 /* ln(maxdouble())*/)))  // for large bt, expm1(bt) overflows
			lnC[ind] = (ln(*pb) :- bt)[ind]
	}
	lnC = lnC :- *plna
	lnlambda = (lnX0 = lnY0 :/ *pgamma) :+ lnC :+ bt
	lnx      = (lnX  = lnY  :/ *pgamma) :+ lnC
}

void asdfEstbernoudiff::transformGradient(real matrix g) {}

// log likelihoods and optionally, 1st & 2nd derivatives thereof
void asdfEstbernoudiff::lnPDF(transmorphic vector params, real colvector lnf, | real scalar todo, real matrix g, struct smatrix h) {
	real scalar N
	real colvector tmp1, tmp2, dlnCdb, dlnfdnu, dlnfdlnlambda, dlnfdlnx, dlngdlna, dlnfdb, nogrowth, lnJacobian, dlngdgamma, 
		d2nulnlambda, d2nulnx, d2lnlambda2, d2lnlambdalnx, d2lnx2, d2lnCdb2, d2nu2, d2lnxdgamma2, gamma2, btdivexpm1bt
	pragma unset todo; pragma unset g; pragma unset dlnfdnu; pragma unset dlnfdlnlambda; pragma unset dlnfdlnx
	pragma unset d2nulnlambda; pragma unset d2nulnx; pragma unset d2lnlambda2; pragma unset d2lnlambdalnx; pragma unset d2lnx2; pragma unset d2nu2

	if (todo == .) todo = 0
	processParams(params)

	if (reflect) {
		if (any(*pnu :< -1)) {
			lnf = J(max((rows(*pnu),rows(lnlambda),rows(lnx))), 1, .)
			if (todo)
				g = J(1, 4, lnf)
      return
		}
	} else if (any(*pnu :> 0)) {
		lnf = J(max((rows(*pnu),rows(lnlambda),rows(lnx))), 1, .)
		if (todo)
			g = J(1, 4, lnf)
		return
	}

	lnf = lnPDFFeller(*pnu, lnlambda, lnx, reflect, todo, dlnfdnu, dlnfdlnlambda, dlnfdlnx, d2nu2, d2nulnlambda, d2nulnx, d2lnlambda2, d2lnlambdalnx, d2lnx2)
	lnJacobian = (lnC :- ln(abs(*pgamma))) :- (1 :- 1 :/ *pgamma) :* lnY  // adjust density for var transforms
	if (Nzeros) lnJacobian[zeroind] = zerofill  // when x = 0, the probability, if non-missing, is a point mass; don't adjust with Jacobian
	lnf = lnf + lnJacobian

	if (todo) {
		if (*pb==0) {
			dlnCdb  = tDelta * (-.5)
      btdivexpm1bt = 0
		} else {
      btdivexpm1bt = editmissing(bt :/ expm1(-bt), 0)        // here, precision is excellent as b->0
			dlnCdb  = (1 :+ btdivexpm1bt) :/ *pb   // here, not so much
			if (N = length(nogrowth = selectindex(abs(bt) :< 1e-14)))
				dlnCdb[nogrowth] = -.5 * (rows(tDelta)>1? tDelta[nogrowth] : J(N, 1, tDelta))
		}

		dlngdlna = dlnfdlnx :+ 1; if (Nzeros) dlngdlna[zeroind] = zerofill
		dlngdlna = dlngdlna  :+ dlnfdlnlambda
		dlnfdb = dlnfdlnlambda :* tDelta + dlnCdb :* dlngdlna
		dlnlambdadgamma = lnX0  :/ -*pgamma
		dlnxdgamma      = lnX   :/ -*pgamma
		dlngdgamma = dlnfdlnlambda :* dlnlambdadgamma + dlnfdlnx :* dlnxdgamma
		tmp1 = (1 :+ lnX) :/ -*pgamma  // derivative of log Jacobian of transform
		if (Nzeros) tmp1[zeroind] = zerofill
		dlngdgamma = dlngdgamma + tmp1
		g = -dlngdlna, dlnfdb, dlnfdnu, dlngdgamma
		transformGradient(g)  // needed only for "bernoudiff2" reparameterization
 
		if (todo > 1) {
			d2lnCdb2 =  (tDelta :* btdivexpm1bt :* bdivexpm1bt :+ 1):/(-*pb :* *pb)
			if (N = length(nogrowth = selectindex(abs(bt) :< .45e-3)))
  			d2lnCdb2[nogrowth] = tDelta2[nogrowth] * -1.5555555555555X-004  // -t^2/12
			d2lnxdgamma2 = (lnX + lnX) :/ (gamma2 = *pgamma * *pgamma)

      h = smatrix(4,4)
			h[1,1].M = (tmp1 = d2lnlambda2 :+ d2lnlambdalnx) + (tmp2 = d2lnlambdalnx :+ d2lnx2)
			h[1,2].M = -dlnCdb:*h[1,1].M :- tDelta:*tmp1
			h[1,3].M = -d2nulnlambda :- d2nulnx
			h[1,4].M =  -dlnlambdadgamma:*tmp1 :- dlnxdgamma:*tmp2
			h[2,2].M =  tDelta2:*d2lnlambda2 :+ ((tDelta+tDelta):*tmp1 :+ dlnCdb:*h[1,1].M):*dlnCdb :+ d2lnCdb2:*dlngdlna
			h[2,3].M =  (dlnCdb:+tDelta):*d2nulnlambda :+ dlnCdb:*d2nulnx
			h[2,4].M =  tDelta:*(dlnlambdadgamma:*d2lnlambda2 :+ dlnxdgamma:*d2lnlambdalnx) :- h[1,4].M :* dlnCdb
			h[3,3].M =  d2nu2
			h[3,4].M =  dlnlambdadgamma:*d2nulnlambda :+ dlnxdgamma:*d2nulnx
			tmp1 = 1 :/ gamma2  :+ d2lnxdgamma2; if (Nzeros) tmp1[zeroind] = zerofill
			h[4,4].M =  dlnlambdadgamma:*dlnlambdadgamma:*d2lnlambda2 :+ 
			           (dlnlambdadgamma+dlnlambdadgamma):*dlnxdgamma:*d2lnlambdalnx :+ 
								  dlnxdgamma:*dlnxdgamma:*d2lnx2 :+ 
									dlnfdlnlambda:*(lnX0+lnX0):/gamma2  :+ dlnfdlnx:*d2lnxdgamma2 :+ tmp1
		}
	}
}

// CDF of observations, based on crude power series implementation
real colvector asdfEstbernoudiff::CDF() {
	real colvector p
  p = CDFFeller(*pnu, lnlambda, lnx, reflect)
	if (*pgamma < 0)
		p = 1 :- p
  return(p)
}

// inverse CDF of quantiles
real colvector asdfEstbernoudiff::invCDF(real colvector p)
	return ((invCDFFeller((*pgamma>0? p : 1:-p), *pnu, lnlambda, reflect) :- lnC) :* *pgamma)


// class for estimation of Feller diffusion as special case of Bernoulli diffusion with gamma=1

class asdfEstfeller extends asdfEstbernoudiff {
	virtual void getScalarParams(), getSmatrixParams(), lnPDF()
	virtual string rowvector getParamEstNames()
	void new()
}

void asdfEstfeller::new()
	params = smatrix(3,1)  // pre-allocated parameter vessel
	
string rowvector asdfEstfeller::getParamEstNames()
	return (("lna", "b", "nu"))

void asdfEstfeller::getScalarParams(real vector params) {
	plna = &(params[1]); pb = &(params[2]); pnu = &(params[3]); pgamma = &1
}

void asdfEstfeller::getSmatrixParams(struct smatrix vector params) {
	plna = &(params[1].M); pb = &(params[2].M); pnu = &(params[3].M); pgamma = &1
}

// log likelihoods and optionally, 1st & 2nd derivatives thereof; based on asdfEstbernoudiff(), taking gamma=1 and dropping it as a parameter
void asdfEstfeller::lnPDF(transmorphic vector params, real colvector lnf, | real scalar todo, real matrix g, struct smatrix h) {
	real scalar N
	real colvector tmp1, dlnCdb, dlnfdnu, dlnfdlnlambda, dlnfdlnx, dlngdlna, dlnfdb, nogrowth, 
		d2nulnlambda, d2nulnx, d2lnlambda2, d2lnlambdalnx, d2lnx2, d2lnCdb2, btdivexpm1bt
	pragma unset todo; pragma unset g; pragma unset dlnfdnu; pragma unset dlnfdlnlambda; pragma unset dlnfdlnx
	pragma unset d2nulnlambda; pragma unset d2nulnx; pragma unset d2lnlambda2; pragma unset d2lnlambdalnx; pragma unset d2lnx2

	if (todo == .) todo = 0
	processParams(params)

	if (reflect) {
		if (any(*pnu :< -1)) {
			lnf = J(max((rows(*pnu),rows(lnlambda),rows(lnx))), 1, .)
			if (todo)
				g = J(1, 3, lnf)
      return
		}
	} else if (any(*pnu :> 0)) {
		lnf = J(max((rows(*pnu),rows(lnlambda),rows(lnx))), 1, .)
		if (todo)
			g = J(1, 3, lnf)
		return
	}

  h = smatrix(3,3); h[3,3].M=J(0,1,0)
	lnf = lnPDFFeller(*pnu, lnlambda, lnx, reflect, todo, dlnfdnu, dlnfdlnlambda, dlnfdlnx, h[3,3].M, d2nulnlambda, d2nulnx, d2lnlambda2, d2lnlambdalnx, d2lnx2)
	if (Nzeros) lnC[zeroind] = zerofill  // when x = 0, the probability, if non-missing, is a point mass; don't adjust with Jacobian
	lnf = lnf + lnC

	if (todo) {
		if (*pb==0) {
			dlnCdb  = tDelta * (-.5)
      btdivexpm1bt = 0
		} else {
      btdivexpm1bt = editmissing(bt :/ expm1(-bt), 0)        // here, precision is excellent as b->0
			dlnCdb  = (1 :+ btdivexpm1bt) :/ *pb   // here, not so much
			if (N = length(nogrowth = selectindex(abs(bt) :< 1e-14)))
				dlnCdb[nogrowth] = -.5 * (rows(tDelta)>1? tDelta[nogrowth] : J(N, 1, tDelta))
		}

		dlngdlna = dlnfdlnx :+ 1; if (Nzeros) dlngdlna[zeroind] = zerofill
		dlngdlna = dlngdlna  :+ dlnfdlnlambda
		dlnfdb = dlnfdlnlambda :* tDelta + dlnCdb :* dlngdlna
		g = -dlngdlna, dlnfdb, dlnfdnu
 
		if (todo > 1) {
			d2lnCdb2 =  (tDelta :* btdivexpm1bt :* bdivexpm1bt :+ 1):/(-*pb :* *pb)
			if (N = length(nogrowth = selectindex(abs(bt) :< .45e-3)))
  			d2lnCdb2[nogrowth] = tDelta2[nogrowth] * -1.5555555555555X-004  // -t^2/12

			h[1,1].M = (tmp1 = d2lnlambda2 :+ d2lnlambdalnx) + (d2lnlambdalnx :+ d2lnx2)
			h[1,2].M = -dlnCdb:*h[1,1].M :- tDelta:*tmp1
			h[1,3].M = -d2nulnlambda :- d2nulnx
			h[2,2].M =  tDelta2:*d2lnlambda2 :+ ((tDelta+tDelta):*tmp1 :+ dlnCdb:*h[1,1].M):*dlnCdb :+ d2lnCdb2:*dlngdlna
			h[2,3].M =  (dlnCdb:+tDelta):*d2nulnlambda :+ dlnCdb:*d2nulnx
			// h[3,3].M =  d2nu2
		}
	}
}

///
/// Sticky Feller routines
///

// Class for computing complex-argument digamma

class clsdigamma {
	real rowvector B
  real colvector k, zero2seven, one2seven
	void new()
	complex matrix psi()
}

void clsdigamma::new() {
  // coefficients and exponents in asymptotic power series, Zhang & Jin (1996), (3.3.14)
  B = -0.5, -1.5555555555555X-004, 1.1111111111111X-007, -1.0410410410410X-008, 1.1111111111111X-008, -1.f07c1f07c1f08X-008, 1.5995995995996X-006, -1.5555555555555X-004, 1.c5e5e5e5e5e5eX-002
  k = -1 \ -2 * (1::8)
  zero2seven = 0::7
  one2seven = 1::7
}

// digamma function for complex argument, based on Zhang and Jin (1996), section 3.3.4
// for simplicity, it *always* applies the recurrence relation 8 times to scaffold to a large (>8) argument appropriate for the asymptotic formula
// this makes it behave smoothly, except when Re(z)=0 toggles usage of the reflection formula
// *** will return large values for non-positive integers but not infinity (missing)
complex matrix clsdigamma::psi(complex matrix z) {
  real scalar i,j; complex scalar lnz; complex matrix retval; real matrix Rez

  Rez = Re(z)
  retval = J(i=rows(z),cols(z),C(.))
  for (;i;i--)
    for (j=cols(z);j;j--)
      if (Rez[i,j] >= 0) {
        lnz = ln(8 + z[i,j])
        retval[i,j] = lnz + B * exp(lnz * k) - quadcolsum(1 :/ (zero2seven :+ z[i,j]))  // asymptotic power series + recurrence relation
      } else {  // exploit reflection formula in case getting near the negative real axis
        lnz = ln(8 - z[i,j])
        retval[i,j] = lnz + B * exp(lnz * k) - quadcolsum(1 :/ (one2seven  :- z[i,j])) - 1.921fb54442d18X+001 /*pi*/ / tan(1.921fb54442d18X+001 /*pi*/ * z[i,j])  // asymptotic power series + recurrence relation + reflection formula
      }
  return(retval)
}


// Class for large-argument Bessel K asymptotic expansion

class clsBesselKAsym {
	real colvector k, twokm1sq
	void new()
	complex rowvector lnK()
}

void clsBesselKAsym::new() {
	k = 0::$BesselKAsym_M
	twokm1sq = 2*k :- 1; twokm1sq = twokm1sq :* twokm1sq
	k = ln(8 * k)
}

// Compute log (K_nu(z)/z^nu) with asymptotic formula for large argument (here, -1<nu<1 ASSUMED) and complex argument z, dlmf.nist.gov/10.40.E2
// doesn't handle z = 0 
// max series length index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// This code only takes first 8 terms of series, which they recommend for abs(z) > 50.
complex rowvector clsBesselKAsym::lnK(real scalar nu, complex rowvector z, complex rowvector lnz, real scalar todo, | complex rowvector dlnKdnu, complex rowvector dlnKdlnz) {
  real scalar twonu, i; complex matrix terms, dtermsdnu, dtermsdlnz; complex rowvector retval; complex colvector fourkmfourksq; complex scalar  _dtermsdnu,  _dtermsdlnz

  if (nu == -.5 | nu == .5)
		if (todo) {
      fourkmfourksq = C(1 :- twokm1sq)
      if (nu == .5) {
        (terms      = J(1, i=cols(z),    ln(fourkmfourksq) - k) :- lnz)[1,] = retval = -lnz - z  // log of first term, including outer multipliers
        (dtermsdnu  = J(1, i        ,  4 :/ fourkmfourksq)         )[1,] = -lnz
        dlnKdlnz = -1 :- z
      } else {
        (terms      = J(1, i=cols(z),    ln(fourkmfourksq) - k) :- lnz     )[1,] = retval = -z  // log of first term, including outer multipliers
        (dtermsdnu  = J(1, i        , -4 :/ fourkmfourksq)                 )[1,] = -lnz
        dlnKdlnz = -z
      }
      dlnKdnu = J(1,i,C(.))
      for (;i;i--) {
        dlnKdnu[i] = quadcolsum(exp(quadrunningsum(terms[,i] :- retval[i], 1) + ln(quadrunningsum(dtermsdnu[,i], 1))))
        retval[i] = 1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ + retval[i]
      }
      return(retval)
    } else
      return(1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :- ((nu + .5) * lnz + z))

	twonu = nu + nu
	(terms = J(1, i=cols(z), ln(C(twonu*twonu :- twokm1sq)) - k) :- lnz)[1,] = (-.5 - nu) * lnz - z  // log of first term, including outer multipliers
  if (todo) {
    (dtermsdnu  = J(1, i, C((4*twonu) :/ (twonu*twonu :- twokm1sq))))[1,] = -lnz
    (dtermsdlnz = J(`=1+$BesselKAsym_M', i, C(-1)))[1,] = (-.5 - nu) :- z
  }

	retval = J(1,i,C(.))
  if (todo) {
    dlnKdnu = dlnKdlnz = J(1,i=cols(z),C(.))
    for (;i;i--) {
      retval[i] = asdfLogSumExpc(quadrunningsum(terms[,i], 1), _dtermsdnu=quadrunningsum(dtermsdnu [,i], 1), _dtermsdlnz=quadrunningsum(dtermsdlnz[,i], 1))
      dlnKdnu [i] = _dtermsdnu
      dlnKdlnz[i] = _dtermsdlnz
      retval[i] = 1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ + retval[i]
    } 
  } else
    for (;i;i--)
      retval[i] = 1.ce6bb25aa1312X-003 + asdfLogSumExpc(quadrunningsum(terms[,i], 1))
	return(retval)
}


// Class for standard power series expansion of Tricomi U function
// Returns log. Normalized so U(a,b,0)=1

class clsUPower {  
  pointer(complex rowvector) scalar palpha
	real scalar nu, paramDirty, lnz, i
  real colvector m, lnm, mlnz, denom1, denom2, mmnu, mpnu, ddenom1dnu, ddenom2dnu
  complex matrix coefs1, coefs2, dcoefs1dalpha, dcoefs2dalpha, dcoefs1dnu, dcoefs2dnu

  void new(), setalpha(), setnu(), setz()
  complex rowvector lnU()
}
	
void clsUPower::new() {
	lnm  = ln(m = 0::$UPower_M)
  dcoefs1dalpha = dcoefs2dalpha = dcoefs1dnu = dcoefs2dnu = coefs1 = coefs2 = J(`=$UPower_M+1', $StickyFeller_N, C(.))
}

void clsUPower::setalpha(complex rowvector _alpha) {
	if (palpha==NULL)
  	palpha = &_alpha
  else if (cols(*palpha) != cols(_alpha) | (*palpha)[1] != _alpha[1])
  	palpha = &_alpha
  paramDirty = 1
}

void clsUPower::setnu(real scalar _nu) {
  if (nu != _nu) {
    nu = _nu
    paramDirty = 1
    denom1 = quadrunningsum(ln(mpnu =  nu :: $UPower_M.5 + nu) + lnm)  // first entry will be 0 since lnm[1]=.
    denom2 = quadrunningsum(ln(mmnu = -nu :: $UPower_M.5 - nu) + lnm)
    
    (ddenom1dnu =   1  :/ mpnu)[1] = 0; ddenom1dnu = quadrunningsum(ddenom1dnu)
    (ddenom2dnu = (-1) :/ mmnu)[1] = 0; ddenom2dnu = quadrunningsum(ddenom2dnu)
  }
}

void clsUPower::setz(real scalar _z, real scalar _lnz) {
	pragma unused _z
  if (lnz != _lnz) {
    lnz = _lnz
    mlnz = m * lnz
  }
}

// return value will be StickyFeller_N long; only firt maxi entries meaningful
complex rowvector clsUPower::lnU(real scalar maxi, complex rowvector lnC, real scalar todo, | complex rowvector dlnCdalpha, complex rowvector dlnCdnu, complex rowvector dlnUdalpha, complex rowvector dlnUdnu, complex rowvector dlnUdlnz) {
  complex scalar _alpha; real scalar _maxi; complex rowvector retval; complex matrix terms; complex colvector dterms2dnu, dterms1dalpha, dterms2dalpha, terms1, terms2, alpham1pm, alpham1mnupm

  _maxi = maxi > $StickyFeller_N? $StickyFeller_N : maxi
	if (paramDirty) i = 1
  for (; i<=_maxi; i++) {
    _alpha = (*palpha)[i]
    alpham1pm    = (_alpha - 1     ) :+ m
    alpham1mnupm = (_alpha - 1 - nu) :+ m

    (terms1 = ln(alpham1pm   ))[1] = 0       // for making log Pochhammer symbols (alpha)_m and (alpha+1-beta)_m; edit m=0 entry to ln 1
    (terms2 = ln(alpham1mnupm))[1] = lnC[i]  // edit m=0 entry to multiplier on second series

		coefs1[,i] = quadrunningsum(terms1) - denom1
		coefs2[,i] = quadrunningsum(terms2) - denom2

    if (todo) {
      (dterms1dalpha =   1  :/ alpham1pm   )[1] = 0
      (dterms2dalpha =   1  :/ alpham1mnupm)[1] = dlnCdalpha[i]
      (dterms2dnu    = (-1) :/ alpham1mnupm)[1] = dlnCdnu   [i]

      dcoefs1dalpha[,i] = quadrunningsum(dterms1dalpha)
      dcoefs2dalpha[,i] = quadrunningsum(dterms2dalpha)
      dcoefs1dnu   [,i] =                               - ddenom1dnu
      dcoefs2dnu   [,i] = quadrunningsum(dterms2dnu   ) - ddenom2dnu
    }
  }

  terms = coefs1 :+ mlnz \ coefs2 :+ (1.921fb54442d18X+001i /*ln -1*/ :+ mmnu * lnz)

  if (todo) {
    dlnUdalpha =   dcoefs1dalpha \ dcoefs2dalpha         
    dlnUdnu    =   dcoefs1dnu    \ dcoefs2dnu :- lnz
    dlnUdlnz   = C(m             \ mmnu)
    retval = asdfLogSumExpc(terms, dlnUdalpha, dlnUdnu, dlnUdlnz)
  } else
    retval = asdfLogSumExpc(terms)
  paramDirty = 0
  return(editmissing(retval,0))  // return power series sum in logs; lnU(0) = 0
}


// Class for asymptotic large-argument expansion of Tricomi U function.
// Returns log. Normalized so U(a,b,0)=1

class clsUAsymArg {
	real scalar alphaDirty, nuDirty, nu, z, lnz, lngammanegnu, dlngammanegnudnu
	real rowvector qa, qb, qc, qd
	complex rowvector lngammaalphamnu
	pointer(complex rowvector) scalar palpha
  class clsdigamma scalar Sdigamma

  void setalpha(), setnu(), setz()
  complex rowvector lnU()
	real rowvector QuarticRoot()
}

// *** Doesn't make a local copy of alpha, just takes by reference
void clsUAsymArg::setalpha(complex rowvector _alpha) {  // maybe can be dispensed with through "extends"
	palpha = &_alpha
	alphaDirty = 1
}

void clsUAsymArg::setnu(real scalar _nu) {
  if (nu != _nu) {
    nu = _nu
		lngammanegnu = lngamma(-nu)
    dlngammanegnudnu = -digamma(-nu)
    nuDirty = 1
  }
}

void clsUAsymArg::setz(real scalar _z, real scalar _lnz) {
	if (lnz != _lnz) {
		  z  =   _z
    lnz  = _lnz
	}
}

// Find highest real root of a quartic equation using Euler's method, mathforum.org/dr.math/faq/faq.cubic.equations.html
real rowvector clsUAsymArg::QuarticRoot(real rowvector a, real rowvector b, real rowvector c, real rowvector d) {
  real rowvector e, f, fourg, h, i, j, a2, h2, h3; complex rowvector alpha, z1, z2, p, q, r; complex matrix x

	    e = b - .375*(a2=a:*a)  // substitute with y = x - b/4
	    f = c + (.125*a2 - .5*b):*a
	fourg = 4*d - a2:*(.046875*a2 - .25*b) - a:*c

	h = .5*e  // auxilliary cubic equation
	i = .0625*(e:*e-fourg)
	j = -.015625*f:*f

	p = (i-(h2=h:*h)/3) / 3  // substite with z = y - h/3
	q = (((2/27)*h2-i/3):*h+j) * .5

	alpha = (sqrt(C(q:*q + p:*p:*p)) - q) :^ (1/3)
	z1 = alpha - p :/ alpha  // roots of p,q equation
	alpha = alpha * (-.5+1.bb67ae8584cabX-001i /*exp(2i*pi()/3)*/)
	z2 = alpha - p :/ alpha

	h3 = h/3
	p = sqrt(z1 - h3)  // square roots of roots of auxilliary cubic
	q = sqrt(z2 - h3)
	r = f:/(p:*q) * -.125

	x = p+q+r\p-q-r\-p+q-r\-p-q+r
	return (colmax(Re(x) :* abs(Im(x):<1e-3))-0.25*a)
}


// log Tricomi U for large argument, using 2F0 expansion, dlmf.nist.gov/13.7.E3
// fails for entries of whose real part is a negative integer because then lngamma(a) wrongly returns missing (bug reported to Stata 9/29/20)
complex rowvector clsUAsymArg::lnU(real scalar maxi, real scalar todo, | complex rowvector dlnUdalpha, complex rowvector dlnUdnu, complex rowvector dlnUdlnz) {
	real rowvector Ream1, ta, ta2, ta3, absa2, Ima2, M; complex rowvector retval; real colvector m, numer; complex colvector dtermsdalpha0, dtermsdnu0, psi, terms, dnumerdnu, dtermsdnu, dtermsdalpha, dtermsdlnz, mpalpham1, mpalpham1mnu; real scalar i, _maxi; complex scalar _alpha

	if (alphaDirty) {
		Ream1 = Re(*palpha) :- 1; absa2 = abs(*palpha); absa2 = absa2:*absa2; Ima2 = Im(*palpha); Ima2 = Ima2:*Ima2; ta3 = (absa2-Ream1-Ream1) :- 1
	}
	if (alphaDirty | nuDirty) {
		lngammaalphamnu = lngamma(*palpha :- nu)

		ta = Ream1 :- (2+nu); ta2 = ta:*ta + Ima2
		qa = 2*(Ream1+ta); qb = ta3+ta2+4*Ream1:*ta; qc = 2*(Ream1:*ta2+ta:*ta3); qd = ta3:*ta2  // coefficients of quartic equation governing when norm of ratio of successive terms is 1
	}
  _maxi = maxi > $StickyFeller_N? $StickyFeller_N : maxi
	M = QuarticRoot(qa[|.\_maxi|], qb[|.\_maxi|] :- z*z, qc[|.\_maxi|], qd[|.\_maxi|])  // indexes about where terms start to permanently increase

	retval = J(1, cols(M), C(.))
  if (todo) {
  	dlnUdnu = dlnUdalpha = dlnUdlnz = retval
    psi = Sdigamma.psi(*palpha :- nu)
    dtermsdnu0 = -psi :- dlngammanegnudnu
    dtermsdalpha0 = psi :- lnz
  }

  if (all(M:<1e6))  // check if too many terms needed for this approach; if so return missing; probably need Temme's expansion for large alpha and z, which I have not implemented
    for (i=cols(M);i;i--) {
    	_alpha = (*palpha)[i]
      m = 0 :: M[i]
      mpalpham1    = m :+ (_alpha - 1)
      mpalpham1mnu = m :+ (_alpha - 1 - nu)
      numer = ln(mpalpham1) + ln(mpalpham1mnu) - ln(m)  // scope for saving and reusing these when only z changes
      (terms = numer :- (1.921fb54442d18X+001i /*pi i*/ + lnz))[1] = lngammaalphamnu[i] - lngammanegnu - _alpha * lnz

      if (todo) {
      	dnumerdnu     = 1 :/ mpalpham1mnu
        (dtermsdnu    =                 -dnumerdnu)[1] = dtermsdnu0
        (dtermsdalpha = 1 :/ mpalpham1 + dnumerdnu)[1] = dtermsdalpha0
        dtermsdlnz    = -_alpha :- m

        retval[i] = asdfLogSumExpc(quadrunningsum(terms), dtermsdnu, dtermsdalpha, dtermsdlnz)

        dlnUdnu   [i] = dtermsdnu
        dlnUdalpha[i] = dtermsdalpha
        dlnUdlnz  [i] = dtermsdlnz
      } else
        retval[i] = asdfLogSumExpc(quadrunningsum(terms))
    }
	alphaDirty = nuDirty = 0
	return (retval)
}

// Class for Abad and Sesma (1997) large-order (large alpha) expansion of (log) Tricomi U function
// Normalized so U(a,b,0)=1

class clsUAbadSesma {
  real scalar nu, z, lnz, nuDirty, beta
	real colvector B, f, dfdnu, lns
  complex colvector g, dgdlnz, Bdiv4i
  pointer(complex rowvector) scalar palpha
  pointer(complex colvector) colvector pcombB
  pointer(real colvector) colvector pcomb, pj, psp1m2j
  class clsBesselKAsym scalar S

  void new(), setalpha(), setnu(), setz()
	complex rowvector lnU()
}

void clsUAbadSesma::new() {
	real scalar s; real colvector j

  // 4^k/k*|B_2k|, B_k Bernoulli numbers; starting with k=1
  B = 1.5555555555555X-001\1.1111111111111X-002\1.0410410410410X-001\1.1111111111111X+001\1.f07c1f07c1f08X+003\1.5995995995995X+007\1.5555555555556X+00b\1.c5e5e5e5e5e5eX+00f\1.86e7f9b9fe6e8X+014\1.a74ca514ca515X+019\1.1975cc0ed7304X+01f\1.c2f0566566567X+024\1.ac572aaaaaaabX+02a\1.dc0b1a5cfbe17X+030\1.31fad7cbf3c00X+037\1.c280563b8bcbdX+03d\1.7892edfdf5556X+044\1.62b8b44651d09X+04b\1.76024c215d22bX+052\1.b6c0dfed2955bX+059\1.1cca39b77b027X+061\1.97212d8cc1040X+068\1.3f0cb06b17e29X+070\1.1101d96823ee1X+078\1.fc474bdd53c20X+07f\1.007db56db95dfX+088\1.17c6dd28a9378X+090\1.48df88a383ad8X+098\1.9f7b3fa37f314X+0a0\1.195c16c40d563X+0a9\1.97922eafb5d17X+0b1\1.3b0a43def5904X+0ba\1.035a171273533X+0c3\1.c5e89f1fd242eX+0cb\1.a575dd47b788aX+0d4\1.9e84a01ae153bX+0dd\1.af26959307253X+0e6\1.d9890bc68eeddX+0ef\1.1231831c6ccbdX+0f9\1.4e5b1c79769adX+102\1.acc2917790916X+10b\1.20bd5e935fc86X+115\1.97f9f572293ebX+11e\1.2e083516087a9X+128\1.d41e650690e84X+131\1.7b59a0f58ffb8X+13b\1.412709736939fX+145\1.1bc48bbc5cc50X+14f\1.0575ec3ab2234X+159\1.f5fe45963cac3X+162\1.f5abf07320dbaX+16c\1.04c0734d5322dX+177\1.19bdb7ca080b5X+181\1.3c2f25da2f119X+18b\1.704a6714a172bX+195\1.bcf2002532ee2X+19f\1.169855a5defafX+1aa\1.6963e57ed41dfX+1b4\1.e54fbbb872117X+1be\1.5125bc119ad89X+1c9\1.e4616e3d8d09aX+1d3\1.679a15f64d5a9X+1de\1.13c15fcb28309X+1e9\1.b49d69397c7aaX+1f3\1.64abed86e0175X+1fe\1.2c8197ead1271X+209\1.05015c2200ecfX+214\1.d32f2c47ba10cX+21e\1.aea4c26df15f6X+229\1.98adca0fa6ca5X+234\1.8f1dcfb3cd633X+23f\1.90f42af6fe074X+24a\1.9e2f8a971c6f1X+255\1.b7c7b687b5bb6X+260\1.dfcb00aa4c2cfX+26b\1.0cd20cf8a2e4dX+277\1.354a937fa88e0X+282\1.6d400c144ce9eX+28d\1.ba92b59047eb7X+298\1.13079715a5f66X+2a4\1.5e8194cdc3f82X+2af\1.c9e5d8fd49582X+2ba\1.32811f8ce591dX+2c6\1.a45eae570d41eX+2d1\1.273c29596492dX+2dd\1.a89a664a7f630X+2e8\1.3888d893937f4X+2f4\1.d6d21fe9b4fb0X+2ff\1.6acfec47defa2X+30b\1.1df41024fe6a3X+317\1.cce916a2f02f1X+322\1.7bbbe1e92e94bX+32e\1.3fbfc0ec34e1dX+33a\1.131bca925b254X+346\1.e39bd06647752X+351\1.b221ac47c86e0X+35d\1.8df13eed18d1cX+369\1.746394821e395X+375\1.63ae13ac3b025X+381\1.5aabe04f8c67bX+38d
  Bdiv4i = B / 4i

  f = J($UAbadSesma_N,1,1)
  g = J($UAbadSesma_N,1,C(1))
  dfdnu = J($UAbadSesma_N,1,0)
  dgdlnz = J($UAbadSesma_N,1,C(0))
  lns = ln(1::$UAbadSesma_N)
  
  pcomb  = pj = J($UAbadSesma_N, 1, &1)
  psp1m2j = J($UAbadSesma_N, 1, &2)
  pcombB = J($UAbadSesma_N, 1, NULL)
  for (s=2;s<$UAbadSesma_N;s++) {
    pcombB[s] = &(*pcomb[s-1] :* Bdiv4i[*pj[s-1] :+ 1])
    j = 0 :: s * .5  // equivalent to j = 0::floor(s/2)
    pj[s] = &(1 :+ j)
    j = j + j
    psp1m2j[s] = &((s + 1) :- j)
    pcomb  [s] = &comb(s, j)
  }
}

void clsUAbadSesma::setalpha(complex rowvector _alpha) {
	if (palpha==NULL)
  	palpha = &_alpha
  else if (cols(*palpha) != cols(_alpha) | (*palpha)[1] != _alpha[1])
  	palpha = &_alpha
}

void clsUAbadSesma::setnu(real scalar _nu) {
  if (nu != _nu) {
    nu = _nu
    nuDirty = 1
  }
}

void clsUAbadSesma::setz(real scalar _z, real scalar _lnz) {
  if (lnz != _lnz) {
    lnz = _lnz
    z   = _z
  }
}

complex rowvector clsUAbadSesma::lnU(real scalar mini, real scalar todo, | complex rowvector dlnUdalpha, complex rowvector dlnUdnu, complex rowvector dlnUdlnz) {
  real colvector j
  complex scalar Karg, Karg2, lnKarg, lnKarg2, lnnegiz, _term, _combfg
  real scalar s, _cross, _comb, onemhalfbeta, twonupsm1
  complex matrix terms
  complex rowvector denom, dlnKdnuL1, dlnKdlnargL1, dlnKdnuL2, dlnKdlnargL2, retval
  complex colvector _g, dlnargdalpha, dlnargdnu, _combf
  pointer (complex rowvector) scalar pK, pKL1, pKL2, pdlnKdnu, pdlnKdlnarg, pdlnKdnuL1, pdlnKdlnargL1, pdlnKdnuL2, pdlnKdlnargL2, pKdlnKdnuL1, pKdlnKdnuL2, pKdlnKdlnargL1, pKdlnKdlnargL2
  pragma unset dlnKdnuL1; pragma unset dlnKdlnargL1; pragma unset dlnKdnuL2; pragma unset dlnKdlnargL2

  if (nuDirty) {
  	beta = 1 + nu
    onemhalfbeta = 1 - .5*beta
    for (s=1;s<$UAbadSesma_N;s++) {
      j = 0 :: s-1
      f[s+1] = onemhalfbeta * (_cross = quadcross(_comb = comb(2*s-1, 2*j) :* B[s:-j], f[|.\s|]))
      if (todo)
        dfdnu[s+1] = onemhalfbeta * quadcross(_comb, dfdnu[|.\s|]) -.5 * _cross
    }
    nuDirty = 0
  }

  Karg = sqrt(Karg2 = (4 * (*palpha)[|mini+1\.|] :- 2 * beta) * z); lnKarg = ln(Karg); lnKarg2 = lnKarg + lnKarg
  pKL2 = &S.lnK(nu, Karg, lnKarg, todo, dlnKdnuL2, dlnKdlnargL2)
  (terms = J($UAbadSesma_N, cols(Karg), C(.)))[1,] = *pKL2  // actually ln(K_nu)
  pKL2 = &exp(*pKL2)

  if (todo) {
    pdlnKdnuL2     = &dlnKdnuL2
    pdlnKdlnargL2  = &dlnKdlnargL2
    pKdlnKdnuL2    = &(*pKL2 :* *pdlnKdnuL2   )
    pKdlnKdlnargL2 = &(*pKL2 :* *pdlnKdlnargL2)

  	dlnargdnu = -.5 * (dlnargdalpha = 1 :/ (2 * (*palpha)[|mini+1\.|] :- beta))
  	// dlnargdnlz = .5
    (dlnUdalpha = terms)[1,] =             dlnKdlnargL2 :* dlnargdalpha
    (dlnUdnu    = terms)[1,] = dlnKdnuL2 + dlnKdlnargL2 :* dlnargdnu
    (dlnUdlnz   = terms)[1,] =             dlnKdlnargL2  * .5
  }

  // unroll loop below for s=1, to prep for recurrence relation
  g[2] = z * -1.5555555555555X-003i
  pKL1 = &S.lnK(beta, Karg, lnKarg, todo, dlnKdnuL1, dlnKdlnargL1)
  terms[2,] = (lnz + lnz -1.cab0bfa2a2002X+000+1.921fb54442d18X+001i /*ln(-1/6)*/) :+ *pKL1  // update lnK using recurrence relation; *pkL1 is actually ln(K_nu)
  pKL1 = &exp(*pKL1)
  
  if (todo) {
    pdlnKdnuL1     = &dlnKdnuL1
    pdlnKdlnargL1  = &dlnKdlnargL1
    pKdlnKdnuL1    = &(*pKL1 :* *pdlnKdnuL1   )
    pKdlnKdlnargL1 = &(*pKL1 :* *pdlnKdlnargL1)

    dgdlnz[2]      = g[2]
    dlnUdalpha[2,] =                          dlnKdlnargL1 :* dlnargdalpha
    dlnUdnu   [2,] = dfdnu[1] :+ (dlnKdnuL1 + dlnKdlnargL1 :* dlnargdnu   )
    dlnUdlnz  [2,] = 2        :+              dlnKdlnargL1 :* .5
  }

  _term = lnnegiz = -1.921fb54442d18X+000i /*log -i*/ + lnz
  twonupsm1 = nu + nu
  for (s=2;s<$UAbadSesma_N;s++) {
  	g[s+1] = z * (*pcombB[s] ' g[*psp1m2j[s-1]])
    if (todo)
      dgdlnz[s+1] = g[s+1] + z * (*pcombB[s] ' dgdlnz[*psp1m2j[s-1]])

    pK = &(ln(denom = *pKL2 + (++twonupsm1) * *pKL1) - lnKarg2)  // recurrence relation for K_nu(z)/z^nu; pK is actually ln(K_nu)
    terms[s+1,] = *pK :+ ((_term = _term + lnnegiz - lns[s]) + ln(_combfg = (_combf = *pcomb[s] :* f[*pj[s]]) ' (_g = g[*psp1m2j[s]])))  // compute lnK(nu+s, Karg) using recurrence relation for K()
    pK = &(denom :/ Karg2)  // now pK is K_nu, not ln K_nu

    if (todo) {
      pdlnKdnu    = &((*pKdlnKdnuL2    + twonupsm1 * *pKdlnKdnuL1    + *pKL1 + *pKL1) :/ denom     )  // derivative recurrence relations for K_nu(z)/z^nu
      pdlnKdlnarg = &((*pKdlnKdlnargL2 + twonupsm1 * *pKdlnKdlnargL1                ) :/ denom :- 2)

    	dlnUdalpha[s+1,] =                                                                                      *pdlnKdlnarg :* dlnargdalpha
    	dlnUdnu   [s+1,] = ((    *pcomb[s] :* dfdnu[*pj[s]]) ' _g                 ) :/ _combfg  :+ (*pdlnKdnu + *pdlnKdlnarg :* dlnargdnu   )
      dlnUdlnz  [s+1,] = (s + (_combf                      ' dgdlnz[*psp1m2j[s]]) :/ _combfg) :+              *pdlnKdlnarg :* .5

      pdlnKdnuL2     = pdlnKdnuL1     // forward shifts for recurrence calculation
      pdlnKdnuL1     = pdlnKdnu
      pdlnKdlnargL2  = pdlnKdlnargL1
      pdlnKdlnargL1  = pdlnKdlnarg
      pKdlnKdnuL2    = pKdlnKdnuL1
      pKdlnKdnuL1    = &(*pKL1 :* *pdlnKdnuL1   )
      pKdlnKdlnargL1 = &(*pKL1 :* *pdlnKdlnargL1)
    }
    pKL2 = pKL1  // forward shifts for recurrence calculation
    pKL1 = pK
  }

  if (todo) {
    retval = asdfLogSumExpc(terms, dlnUdalpha, dlnUdnu, dlnUdlnz)
    dlnUdnu  = dlnUdnu  :+ (1.62e42fefa39efX-001 /*ln(2)*/ + digamma(-nu)) 
    dlnUdlnz = dlnUdlnz :+ z*.5                                            
  } else
    retval = asdfLogSumExpc(terms)

  return(retval :+ (beta * 1.62e42fefa39efX-001 /*ln(2)*/ + z*.5 - lngamma(-nu)))
}


// Class for computing sticky Feller PDF

class clsStickyFeller {
	real colvector lnBesselKPower_0toM, lnBesselKPower_2Mto0, ot, ox, oy, zerofill, zeroind, missingfill, missingfillg, lnx, lny, Nuniq, _mu, JN10
	real matrix uniqdata
	pointer (real colvector) scalar px, py, pt
	real scalar N, Nt
	complex rowvector tlambdalnu, u2
	complex colvector phi, dphidalpha, dphidnu, dphidlnz
	complex matrix lambda
  class clsUPower scalar SUPower  // instantiate these objectsion when clsStickyFeller is instantiated--just once per estimation
  class clsUAbadSesma scalar SUAsymParam
  class clsUAsymArg scalar SUAsymArg
	class clsBesselKAsym scalar SBesselKAsym
  class clsdigamma scalar Sdigamma
  class asdfEstfeller scalar SEstFeller  // (should be able to access non-sticky Feller functionality more elegantly with class reorganization)

	real colvector lnStickyFeller()
	complex rowvector lnBesselKPower()  // for consistency, should be given its own class
	void new(), setData()
}

void clsStickyFeller::new() {
  complex rowvector u

  SEstFeller.setFixedParam("reflect", 0)  // prep for absorbing term

  lnBesselKPower_0toM = 0::$BesselKPower_M  // stuff for power series
  lnBesselKPower_2Mto0 = lnBesselKPower_0toM * 2
	
	u = C(1, (.5 .. $StickyFeller_N-.5) * $StickyFeller_h)  // really (1 + ui)
	u2 = u :* u
	tlambdalnu = ln(u) + $StickyFeller_mu * u2 :- 1.16c0cdab51479X+001 // + ln 2*2*h/2pi for Mellin integral; tack on here for efficiency. strofreal() always includes +-
}

void clsStickyFeller::setData(real colvector t, real colvector x, real colvector y) {
	real colvector tID, txyID, o, invo; real matrix data, tinfo, txyinfo
	pragma unset tID; pragma unset txyID

	SEstFeller.setData(t, x, y)

  px = &x; py = &y; pt = &t
	
	zerofill = J(length(zeroind = selectindex(y:==0)), 1, 0)
	missingfill  = J(length(zeroind), 1, .)
	missingfillg = J(length(zeroind), 3, 0)
	lnx = ln(x); lny = ln(y)

	N = rows(t)
  JN10 = J(N,1,0)
	data = t,x,lnx \ t,y,lny
	_collate(data, o = order(data, (1,2)))  // sort observations by t, x/y
	tinfo   = asdfpanelsetup(data,  1   ,   tID)
	txyinfo = asdfpanelsetup(data, (1,2), txyID)
	uniqdata = data[txyinfo[,1],]  // unique t-x/y combinations
	Nuniq = rows(txyinfo)
	Nt = rows(tinfo)  // unique values of t
	invo = invorder(o)
	ot =   tID[invo[|.  \N|]]  // map from original observations to unique t's
	ox = txyID[invo[|.  \N|]]  // map from original observations to unique t-x's
	oy = txyID[invo[|N+1\.|]]  // map from original observations to unique t-y's
	dphidalpha = dphidnu = dphidlnz = phi = J(Nuniq, $StickyFeller_N, C(0))

	_mu = $StickyFeller_mu :/ data[tinfo[,1],1]
/*"$StickyFeller_mu,data[tinfo[,1],1]"
 $StickyFeller_mu,data[tinfo[,1],1]*/
	lambda = _mu * u2   // "mu" as in Garrappa (2015), p. 1364, not stickiness parameter
}


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007), Weideman (2010), parabolic contour method to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real scalar lna, real scalar b, real scalar nu, real scalar lnmu, real scalar todo, | real matrix g) {
  real colvector retval, lnp0, lnm, z, lnz, dlnmdlna, dlnmdb, dlnmdnu, dlnmdlnmu
	real scalar a, hi, ba, i, j, lngammanegnu, beta, p, absb, anypower, anyasymparam, anyasymarg, oldt, lnba, absbdiv_mu, c, floorp, sqBessel
  complex matrix dlna, db, dnu, dlnmu, alpha, alphamnu, lngammaalphamnu, lnK_mu, lnC, dlnCdb, dlnCdnu, dlnCdalpha, K_mu, KC, dlnKdlna, dlnKdb, dlnKdnu, dlnKdlnmu, terms, dphidb
	complex rowvector lnUAsymParam, lnUSmallalpha, alphaj, lnalphaj, lambdaj, lnCj, dlnCdalphaj, dlnCdnuj, dlnUdalphaSmall, dlnUdnuSmall, dlnUdlnzSmall, dlnUdalphaBig, dlnUdnuBig, dlnUdlnzBig
  complex colvector lnsticky
	pragma unset dlnUdalphaSmall; pragma unset dlnUdnuSmall; pragma unset dlnUdlnzSmall; pragma unset dlnUdalphaBig; pragma unset dlnUdnuBig; pragma unset dlnUdlnzBig; pragma unset lnp0
  pragma unset oldt  // will exploit that initial value is missing

	if (nu < -1 | nu > 0)
		return (.)

  SEstFeller.lnPDF(lna\b\nu, lnp0, todo, g)
  if (nu == 0 | nu == -1) {
    if (todo)
      g = g, J(rows(g),1,0)  // add 0's for derivatives w.r.t. lnmu
		return (lnp0)
  }    
 
	if (rows(zerofill)) {
    lnp0[zeroind] = missingfill  // where y=0, "absorbing term" drops out
    if (todo)
      g[zeroind,] = missingfillg
  }

	beta = 1 + nu
	lngammanegnu = lngamma(-nu)
	ba = (absb = abs(b)) / (a = exp(lna))
	_editmissing(lnm = nu * lny - ba * (b < 0? *py : *px) :- lna, -lnmu)  // log speed measure; takes special value of 1/mu when y = 0; b>0 modification is a trick to multiply by exp(-b/a*(x+y))
  if (todo) {
    dlnmdlna  = (lny :< .) :* (ba * (b < 0? *py : *px) :- 1)
    dlnmdb    = (lny :< .) :* (b < 0? *py / a : *px / (-a))
    dlnmdnu   = editmissing(lny, 0)
    dlnmdlnmu = -!*py  // -1 if y is 0; 0 otherwise
  }

	if (sqBessel = absb < 1e-10) {  // b=0 case. Sort of squared Bessel, except that would mean a=2. Would need a different approach if b could vary by observation
		lnC = (lngamma(beta) - ln(-nu) - lngammanegnu) :- nu * (lnba = ln(lambda) :- lna)
		lnK_mu = -ln(u2 / exp(lnmu) :- (nu :/ _mu) :* exp(lnC))   // K in G = K * phi(x) * phi(y); mu is stickiness parameter; _mu is mu in Garrappa (2015), p. 1364 (part of Jacobian of parabolic path function)

    if (todo) {
      dlnCdnu = digamma(beta) + digamma(1 - nu)
      KC = exp(lnC) :* (K_mu = exp(lnK_mu) / _mu) :* nu
    }

		j = Nt + 1
		for (i=Nuniq; i; i--) {
			if (uniqdata[i,1] != oldt) {
				oldt = uniqdata[i,1]
				lambdaj = lambda[--j,]
			}

			if (z = uniqdata[i,2]) {
        alphaj = sqrt(lambdaj :* ((4/a) * z)); lnalphaj = ln(alphaj)

        p = editmissing(round( sqrt((20.25 * a / _mu[j]) / z :- 1) / $StickyFeller_h), 0)  // indexes of parabola points where BesselK argument passes 9, cut-off between power series & asymptotic; 20.25 = 9^2/4 
        phi[i,] = (p < 1?                                                                                                  SBesselKAsym.lnK(nu, alphaj,          lnalphaj         , todo, dlnUdnuSmall, dlnUdalphaSmall)  : 
                          (p >= $StickyFeller_N? lnBesselKPower(nu, lnalphaj       , todo, dlnUdnuSmall, dlnUdalphaSmall)  :
                                                 lnBesselKPower(nu, lnalphaj[|.\p|], todo, dlnUdnuSmall, dlnUdalphaSmall), SBesselKAsym.lnK(nu, alphaj[|p+1\.|], lnalphaj[|p+1\.|], todo, dlnUdnuBig  , dlnUdalphaBig))) :+
                  (beta * 1.62e42fefa39efX-001 /*ln 2*/ - lngammanegnu)

        if (todo) {
          if (p >= 1 & p < $StickyFeller_N) {
            dlnUdnuSmall    = dlnUdnuSmall   , dlnUdnuBig
            dlnUdalphaSmall = dlnUdalphaSmall, dlnUdalphaBig
          }
          dphidnu [i,] = dlnUdnuSmall :+ (1.62e42fefa39efX-001 /*ln 2*/ + digamma(-nu))
          dphidlnz[i,] = .5 * dlnUdalphaSmall
        }
      }
		}
  } else {  // non-squared-Bessel
		SUPower.setnu    (nu)
		SUAsymParam.setnu(nu)
		SUAsymArg.setnu  (nu)
		alpha = lambda / -b; if (b > 0) alpha = beta :- alpha

		lngammaalphamnu = lngamma(alphamnu = alpha :- nu)
		lnC = (lngamma(beta) - ln(-nu) - lngammanegnu) :+ (lngammaalphamnu - lngamma(alpha))
		lnK_mu = -ln(u2 / exp(lnmu) :- (nu :/ _mu) :* exp(lnC :- nu * (lnba = ln(absb)-lna)))  // ln K*_mu in G = K * phi(x) * phi(y); mu is stickiness parameter; _mu is mu as in Garrappa (2015), p. 1364, part of Jacobian of parabolic path function
    if (todo) {
      dlnCdb   = lambda / (-b * absb) :* (dlnCdalpha = Sdigamma.psi(alphamnu) - Sdigamma.psi(alpha))
      dlnCdnu  = (digamma(beta) + digamma(1 - nu)) :- Sdigamma.psi(alphamnu)
      KC = exp(lnC) :* (K_mu = exp(lnK_mu) / _mu) * (ba^-nu * nu)
      dlnKdb = KC :* (dlnCdb  :- nu / b)
    }

		j = Nt + 1
		for (i=Nuniq; i; i--) {
			if (uniqdata[i,1] != oldt) {
				oldt = uniqdata[i,1]
				alphaj = alpha[--j,]
				absbdiv_mu = absb / _mu[j]
				dlnCdalphaj = dlnCdnuj = lnCj = C(.)
				SUPower.setalpha(alphaj)
				SUAsymParam.setalpha(alphaj)
				SUAsymArg.setalpha(alphaj)
			}

			z = ba * uniqdata[i,2]  // z = b/a *x or b/a * y
      if (z) {  // if z==0, ln phi=ln 1, as initialized; ditto for derivatives
        lnz = lnba :+ uniqdata[i,3]

        // indexes in lambda after which switch from power or asymptotic-argument to asymptotic-parameter; min value is 0, based on benchmarking against mpmath's hyperu()
        hi = z > 15.460074  // for low z, use large-argument asymptotic rather than standard power series representation?
        p = editmissing(.5 + sqrt( (hi? (z * `=exp(-1.446009 /  1.102413)') ^  1.102413 : 
                                        (z * `=exp( 3.582745 / -.7823305)') ^ -.7823305   ) * absbdiv_mu - 1) / $StickyFeller_h, 0)  // index of last point on parabola handled with power series representation
        floorp = floor(p)
        anypower = hi==0 & p>0; anyasymarg = hi & floorp; anyasymparam = p < $StickyFeller_N
        if (anypower) {
          SUPower.setz(z, lnz)
          if (lnCj==.) {
          	lnCj=lnC[j,]
            if (todo) {
            	dlnCdalphaj = dlnCdalpha[j,]
              dlnCdnuj    = dlnCdnu   [j,]
            }
          }
          lnUSmallalpha = SUPower.lnU(floorp+1, lnCj, todo, dlnCdalphaj, dlnCdnuj, dlnUdalphaSmall, dlnUdnuSmall, dlnUdlnzSmall)
        } else if (anyasymarg) {
          SUAsymArg.setz(z, lnz)
          lnUSmallalpha = SUAsymArg.lnU(floorp+1, todo, dlnUdalphaSmall, dlnUdnuSmall, dlnUdlnzSmall)
        }

        if (anyasymparam) {
          SUAsymParam.setz(z, lnz)
          lnUAsymParam = SUAsymParam.lnU(floorp, todo, dlnUdalphaBig, dlnUdnuBig, dlnUdlnzBig)
        }

        if (p > 0 & anyasymparam) {  // if using both power series/large-argument and large-parameter, smooth cross-over with a weighted average to reduce discontinuities w.r.t parameters
          c = p - floorp
          lnUAsymParam[1] = (1 - c) * lnUAsymParam[1] + c * lnUSmallalpha[floorp+1]
          if (todo) {
            dlnUdalphaBig[1] = (1 - c) * dlnUdalphaBig[1] + c * dlnUdalphaSmall[floorp+1]
            dlnUdnuBig   [1] = (1 - c) * dlnUdnuBig   [1] + c * dlnUdnuSmall   [floorp+1]
            dlnUdlnzBig  [1] = (1 - c) * dlnUdlnzBig  [1] + c * dlnUdlnzSmall  [floorp+1]
          }
        }

        phi[i,] = floorp? (p >= $StickyFeller_N? lnUSmallalpha : lnUSmallalpha[|.\floorp|], lnUAsymParam) : lnUAsymParam
        if (todo) {
          dphidalpha[i,] = floorp? (p >= $StickyFeller_N? dlnUdalphaSmall : dlnUdalphaSmall[|.\floorp|], dlnUdalphaBig) : dlnUdalphaBig
          dphidnu   [i,] = floorp? (p >= $StickyFeller_N? dlnUdnuSmall    : dlnUdnuSmall   [|.\floorp|], dlnUdnuBig   ) : dlnUdnuBig
          dphidlnz  [i,] = floorp? (p >= $StickyFeller_N? dlnUdlnzSmall   : dlnUdlnzSmall  [|.\floorp|], dlnUdlnzBig  ) : dlnUdlnzBig
        }
      }
		}

    if (todo) {
      dphidb  = dphidalpha :* lambda /(-b * absb) :+ dphidlnz / b  // from partial to total derivatives 
      dphidnu = (b < 0? dphidnu : dphidalpha + dphidnu)
    }
	}

  terms = (tlambdalnu :+ lnK_mu)[ot,] + phi[ox,] + phi[oy,]

  if (todo) {
    dlnKdlna  = KC * nu
    dlnKdnu   = KC :* (dlnCdnu :+ 1 / nu :- lnba)
    dlnKdlnmu = K_mu :* lambda / exp(lnmu)

    dlna  =   dlnKdlna [ot,] - dphidlnz[ox,] - dphidlnz[oy,]  // dphidlna = -dphidlnz
    if (sqBessel==0)
      db    = dlnKdb   [ot,] + dphidb  [ox,] + dphidb  [oy,]
    dnu   =   dlnKdnu  [ot,] + dphidnu [ox,] + dphidnu [oy,]
    dlnmu =   dlnKdlnmu[ot,]

    lnsticky = lnm + asdfLogRowSumExpc(terms, dlna, db, dnu, dlnmu)  // add log speed measure: p = m * f

    dlna  = g[,1],                        dlna  + dlnmdlna
    db    = g[,2], (sqBessel? C(dlnmdb) : db    + dlnmdb   )
    dnu   = g[,3],                        dnu   + dlnmdnu
    dlnmu = JN10 ,                        dlnmu + dlnmdlnmu

    retval = Re(asdfLogRowSumExpc((lnp0, lnsticky), dlna, db, dnu, dlnmu))  // Add asborbing & sticky terms. In theory can take Re() sooner, but in extreme cases imprecision makes Bromwich sum <0, throwing us into C and generating discontinuity
    g = Re(dlna), Re(db), Re(dnu), Re(dlnmu)
    return(retval)
  }
  return(Re(asdfLogRowSumExpc((lnp0, asdfLogRowSumExpc(terms) + lnm))))  // add log speed measure: p = m * f, then add absorbing term
}


// Compute log (K_nu(z) / z^nu) with power series formula and complex argument z
// series length from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
complex rowvector clsStickyFeller::lnBesselKPower(real scalar nu, complex rowvector lnz, real scalar todo, | complex rowvector dnu, complex rowvector dlnz) {
  real colvector InuDenom, InegnuDenom, dInuDenomdnu, dInegnuDenomdnu, t1, t2; complex matrix twomlnhalfz, terms; complex rowvector lnhalfz, retval, t3; real scalar pinu

  (InegnuDenom = ln((t1 = -nu::$BesselKPower_M.5-nu) :* lnBesselKPower_0toM))[1] = lngamma(1-nu)
  (InuDenom    = ln((t2 =  nu::$BesselKPower_M.5+nu) :* lnBesselKPower_0toM))[1] = lngamma(1+nu)

	lnhalfz = lnz :- 1.62e42fefa39efX-001 /*ln .5*/
	twomlnhalfz = lnBesselKPower_2Mto0 * lnhalfz
  terms = ((twomlnhalfz :+ nu * (t3 = 1.62e42fefa39efX-000 /*ln 4*/ :- 2*lnz)) :- quadrunningsum(InegnuDenom)) \ 
	        ( twomlnhalfz                                                        :- quadrunningsum(   InuDenom)) :+ 1.921fb54442d18X+001i /*pi i*/

  pinu = nu * 1.921fb54442d18X+001 /*pi*/
  if (todo) {
    (dInegnuDenomdnu =   1  :/ t1)[1] =  digamma(1-nu)  // negated for efficiency
    (dInuDenomdnu    = (-1) :/ t2)[1] = -digamma(1+nu)
    dnu  = t3 :+ J(1,cols(lnz),quadrunningsum(dInegnuDenomdnu)) \
                 J(1,cols(lnz),quadrunningsum(   dInuDenomdnu))
    dlnz = C(lnBesselKPower_2Mto0 :- 2*nu \ 
	           lnBesselKPower_2Mto0            )
    retval = asdfLogSumExpc(terms, dnu, dlnz) :+ (1.ce6bb25aa1315X-002 /*ln pi/2*/ - nu * 1.62e42fefa39efX-001 /*ln 2*/ -                           ln(C(sin(pinu))))
    dnu    = dnu                              :- (                                        1.62e42fefa39efX-001 /*ln 2*/ +  1.921fb54442d18X+001 /*pi*/ / tan(pinu))
  } else
  	retval = asdfLogSumExpc(terms) :+ (1.ce6bb25aa1315X-002 /*ln pi/2*/ - nu * 1.62e42fefa39efX-001 /*ln 2*/ - ln(C(sin(pinu))))
  return(retval)
}

//
// Class for estimating Bernoulli SDE via s, B, delta, lnsig parameterization
//

class asdfEstbernoudiff2 extends asdfEstbernoudiff {
	pointer (real colvector) scalar plnsig, pdelta, ps, pB

	virtual string rowvector getParamEstNames()
	virtual void getScalarParams(), getSmatrixParams(), transformGradient()
	virtual void transformParams()
}

string rowvector asdfEstbernoudiff2::getParamEstNames()
	return (("s", "B", "delta", "lnsig"))
	
void asdfEstbernoudiff2::getScalarParams(real vector params) {
	ps = &(params[1]); pB = &(params[2]); pdelta = &(params[3]); plnsig = &(params[4])
	transformParams()
}

void asdfEstbernoudiff2::getSmatrixParams(struct smatrix vector params) {
	ps = &(params[1].M); pB = &(params[2].M); pdelta = &(params[3].M); plnsig = &(params[4].M)
	transformParams()
}

void asdfEstbernoudiff2::transformParams() {
	real colvector t
	t = *plnsig + ln(abs(*pB))
	plna = &(t + t  :- 1.62e42fefa39efX-001 /*ln(2)*/)
	pb = &(-*pB :* *pdelta)
	pnu = &((1 :- 2* *ps :/ exp(*plnsig + *plnsig)) :/ *pB)
	pgamma = &(-1 :/ *pB)
}

void asdfEstbernoudiff2::transformGradient(real matrix g) {
	real colvector t
	t = -2 :* exp(-*plnsig - *plnsig)
	g = t :/ *pB :* g[,3], ///
			(2* g[,1] + *pb :* g[,2] - *pnu :* g[,3] - *pgamma :* g[,4]) :/ *pB, ///
	    -*pB :* g[,2], ///
			2*g[,1] - (t+t) :* *ps :/ *pB :*g[,3]
}


//
// Class for estimating geometric Brownian motion
//

class asdfEstgbm extends asdfEst {
	real colvector a
	pointer (real colvector) scalar pb

	virtual void lnPDF(), setFixedParam(), CDF(), getScalarParams(), getSmatrixParams()
	virtual string rowvector getParamEstNames(), getParamFixedNames()
	void new()
}

void asdfEstgbm::new() params = smatrix(2,1)  // pre-allocated parameter vessel

string rowvector asdfEstgbm::getParamEstNames()
	return (("lna", "b"))

string rowvector asdfEstgbm::getParamFixedNames()
	return ("reflect")

void asdfEstgbm::getScalarParams(real vector params) {
	a = exp(params[1]); pb = &(params[2])
}

void asdfEstgbm::getSmatrixParams(struct smatrix vector params) {
	a = exp(params[1].M); pb = &(params[2].M)
}

void asdfEstgbm::setFixedParam(string scalar name, numeric value) {
	if (name == "reflect") reflect = value
}

void asdfEstgbm::lnPDF(transmorphic vector params, real colvector lnf, | real scalar todo, real matrix g, struct smatrix h) {
	real colvector z, at

	if (eltype(params)=="real")
		getScalarParams(params)
	else
		getSmatrixParams(params)
	if (todo == .) todo = 0

	z = lnY :- (lnY0 :+ (*pb :- a) :* tDelta)
	at = a :* tDelta
	lnf = lnnormalden(z, sqrt(at+at)) - lnY

	if (todo) {
		g = .25*z:*z:/at :- .5 * (z:+1) , z :/ (a+a)
		if (todo > 1) {
			h = smatrix(2,2)
			h[1,1].M = -.5 * (at        :- z    :+ .5*z:*z:/at)
			h[1,2].M =  .5 * (tDelta    :+ z:/a               )
			h[2,2].M = -.5 *  tDelta:/a
		}
	}
}


//
// class for estimating sticky Feller
//

class asdfEststickyfeller extends asdfEst {
	pointer (real colvector) scalar pb, plna, pnu, plnmu, pY0, pY, ptDelta
	class clsStickyFeller scalar S

	void new(), setData()
	virtual void lnPDF(), setFixedParam(), getScalarParams(), getSmatrixParams(), processParams()
	virtual string rowvector getParamEstNames(), getParamFixedNames()
	real scalar getlf()
}

real scalar asdfEststickyfeller::getlf() return(1)  // lf1 estimator

// _Y0, _Y, _tDelta = inital & final values, and time gaps
void asdfEststickyfeller::setData(real colvector tDelta, real colvector Y0, real colvector Y) {
	pY0 = &Y0; pY  = &Y; ptDelta = &tDelta
	S.setData(tDelta, Y0, Y)
}

string rowvector asdfEststickyfeller::getParamEstNames()
	return (("lna", "b", "nu", "lnmu"))

string rowvector asdfEststickyfeller::getParamFixedNames()
	return ("")

void asdfEststickyfeller::new()
	params = smatrix(4,1)  // pre-allocated parameter vessel

void asdfEststickyfeller::getScalarParams(real vector params) {
	plna = &(params[1]); pb = &(params[2]); pnu = &(params[3]); plnmu = &(params[4])
}

void asdfEststickyfeller::getSmatrixParams(struct smatrix vector params) {
	plna = &(params[1].M); pb = &(params[2].M); pnu = &(params[3].M); plnmu = &(params[4].M)
}

void asdfEststickyfeller::processParams(transmorphic vector params) {
	if (eltype(params)=="real")
		getScalarParams(params)
	else
		getSmatrixParams(params)
}

void asdfEststickyfeller::lnPDF(transmorphic vector params, real colvector lnf, | real scalar todo, real matrix g, struct smatrix h) {
	pragma unset todo; pragma unset g; pragma unset h
	processParams(params)
	lnf = S.lnStickyFeller(*plna, *pb, *pnu, *plnmu, todo, g)
}


//
// class for sticky squared Bessel, as special case of sticky Feller
//

class asdfEststickysqbessel extends asdfEststickyfeller {
	void new()
	virtual void getScalarParams(), getSmatrixParams(), lnPDF()
	virtual string rowvector getParamEstNames()
}

string rowvector asdfEststickysqbessel::getParamEstNames()
	return (("nu", "lnmu"))

void asdfEststickysqbessel::new() {
	params = smatrix(2,1)  // pre-allocated parameter vessel
	plna = &1.62e42fefa39efX-001 /*ln 2*/
	pb = &0
}

void asdfEststickysqbessel::getScalarParams(real vector params) {
	pnu = &(params[1]); plnmu = &(params[2])
}

void asdfEststickysqbessel::getSmatrixParams(struct smatrix vector params) {
	pnu = &(params[1].M); plnmu = &(params[2].M)
}


void asdfEststickysqbessel::lnPDF(transmorphic vector params, real colvector lnf, | real scalar todo, real matrix g, struct smatrix h) {
	pragma unset todo; pragma unset g; pragma unset h
	processParams(params)
	lnf = S.lnStickyFeller(*plna, *pb, *pnu, *plnmu, todo, g)
  if (todo) g = g[|.,3\.,4|]
}

// lf2 likelihood evaluator for all estimation classes
// 1st argument is a moptimize() obect whose 1st userinfo item is the asdfEst obect
function asdflf2(transmorphic scalar M, real scalar todo, real rowvector b, real colvector lnf, real matrix g, real matrix H) {
	class asdfEst scalar S; struct smatrix matrix h, params; real scalar d, i, j; real matrix t
	pragma unset h
	S = moptimize_util_userinfo(M, 1)
	d = rows(params = S.getParams())
	for (i=d;i;i--)
		params[i].M = moptimize_util_xb(M, b, i)
	S.lnPDF(params, lnf, todo, g, h)

	if (todo > 1 & lnf[1] < .) {
		for (i=d; i; i--) {
			h[i,i].M = moptimize_util_matsum(M, i, i, h[i,i].M, 1)
			for(j=i-1; j; j--)
				h[i,j].M = (h[j,i].M = moptimize_util_matsum(M, j, i, h[j,i].M, 1))'
		}
		for (i=1; i<=d; i++) {
			t = h[i,d].M
			for (j=d-1; j; j--)
				t = h[i,j].M , t
			H = i==1? t : H \ t
		}
  }
}

numeric rowvector asdfLogSumExp(numeric matrix x) {
	real rowvector shift; real colvector limits, minn, maxx; real matrix t
	if (rows(x)==0) return(J(1,cols(x),0))
	if (rows(x)==1) return(x)
	limits = -1.61b2bdd7abcd2X+009 /*ln(smallestdouble()) + 1*/ \ 1.620b76e3a7b61X+009 /*ln(maxdouble()) - 1*/ - rows(x)
	t = limits :- colminmax(Re(x))  // shift just enough to prevent underflow in exp(), but if necessary even less to avoid overflow
	minn = t[1,]; maxx = t[2,]
	t = minn :* (minn :> 0) - maxx
	shift = t :* (t :< 0) + maxx // parallelizes better than rowminmax()
	return (any(shift)? ln(quadcolsum(exp(x :+ shift))) - shift : ln(quadcolsum(exp(x))))
}




// logsumexp specifically for complex arguments, treating missing as log of 0. Accepts up to 4 optional derivative matrices and replaces them with corresponding derivatives of result
complex rowvector asdfLogSumExpc(complex matrix x, | complex matrix d1, complex matrix d2, complex matrix d3, complex matrix d4) {
	real rowvector shift; real colvector limits, minn, maxx; real matrix t; complex matrix expx; complex rowvector expretval, retval; real scalar i, anyshift; pointer (complex matrix) colvector pd
	if (rows(x)==0) return(J(1,cols(x),C(0)))
	if (rows(x)==1) return(x)
	limits = -1.61b2bdd7abcd2X+009 /*ln(smallestdouble()) + 1*/ \ 1.620b76e3a7b61X+009 /*ln(maxdouble()) - 1*/ - rows(x)
	t = limits :- colminmax(Re(x))  // shift just enough to prevent underflow in exp(), but if necessary even less to avoid overflow
	minn = t[1,]; maxx = t[2,]
	t = minn :* (minn :> 0) - maxx
	shift = t :* (t :< 0) + maxx // parallelizes better than rowminmax()
  if (args() > 1) {
    retval = (anyshift=any(shift))? ln(quadcolsum(expx = exp(x :+ shift))) - shift : ln(quadcolsum(expx = exp(x)))
    pd = &d1 \ &d2 \ &d3 \ &d4
    expretval = anyshift? exp(retval + shift) : exp(retval)
    for (i=args()-1;i;i--)
      if (cols(*pd[i]))  // skip empty provided-derivative matrices
        *pd[i] = quadcolsum(expx :* *pd[i]) :/ expretval
  } else
    retval = any(shift)? ln(quadcolsum(exp(x :+ shift))) - shift : ln(quadcolsum(exp(x)))
  return(retval)
}

// same, but rowwise
complex colvector asdfLogRowSumExpc(complex matrix x, | complex matrix d1, complex matrix d2, complex matrix d3, complex matrix d4) {
	real colvector shift; real rowvector limits, minn, maxx; real matrix t; complex matrix expx; complex colvector expretval, retval; real scalar i, anyshift; pointer (complex matrix) colvector pd
	if (cols(x)==1) return(x)
	limits = -1.61b2bdd7abcd2X+009 /*ln(smallestdouble()) + 1*/ , 1.620b76e3a7b61X+009 /*ln(maxdouble()) - 1*/ - cols(x)
	t = limits :- rowminmax(Re(x))
	minn = t[,1]; maxx = t[,2]
	t = minn :* (minn :> 0) - maxx
	shift = t :* (t :< 0) + maxx
  if (args() > 1) {
    retval = (anyshift=any(shift))? ln(quadrowsum(expx = exp(x :+ shift))) - shift : ln(quadrowsum(expx = exp(x)))
    pd = &d1 \ &d2 \ &d3 \ &d4
    expretval = anyshift? exp(retval + shift) : exp(retval)
    for (i=args()-1;i;i--)
      if (rows(*pd[i]))  // skip empty provided-derivative matrices
        *pd[i] = quadrowsum(expx :* *pd[i]) :/ expretval
  } else
    retval = any(shift)? ln(quadrowsum(exp(x :+ shift))) - shift : ln(quadrowsum(exp(x)))
  return(retval)
}


// like panelsetup() but can group on multiple columns, like sort(), and faster. But doesn't take minobs, maxobs arguments.
// Takes third argument, a matrix in which to store standardized ID variable, starting from 1
real matrix asdfpanelsetup(real matrix X, real rowvector cols, real colvector ID) {
	real matrix info; real scalar i, N; real scalar p; real rowvector t, id
	N = rows(X)
	info = J(N, 2, N)
  ID = J(N, 1, 1)
	info[1,1] = p = 1
	id = X[1, cols]
	for (i=2; i<=N; i++) {
		if ((t = X[i,cols]) != id) {
			info[  p,2] = i - 1
			info[++p,1] = i
			id = t
		}
		ID[i] = p
	}
	return (info[|.,.\p,.|])
}


// exponentiate a column, treating missing as log of 0
real colvector asdfExpNegInfty(real colvector lnv) {
	real colvector overflow, retval; real scalar t
	retval = editmissing(exp(lnv), 0)
	overflow = selectindex(lnv :>= 708.98956571282406 :& lnv :< .)
	if (t=length(overflow))
		retval[overflow] = J(t, 1, .)
	return (retval)
}


// compute log of f^*_(+/-) (x; lambda, nu) and optional 1st and second derivatives with respect to inputs
// lambda and x are to be passed in logs for computational convenience but the derivatives are still w.r.t. x and lambda
// lnlambdav or lnxv = . interpreted as ln 0
// when x = 0 and there is a mass point at x=0, returns log of that
// reflect == 0 for absorbing diffusion; reflecting otherwise
// Depending on inputs, uses ascending power series and large-argument and large-order asymptotic series for modified Bessel function of first kind. (dlmf.nist.gov, eqs 10.25.2, 10.40.1, 10.41.1)
real colvector lnPDFFeller(real colvector _nuv, real colvector lnlambdav, real colvector lnxv, | real scalar reflect, real scalar _todo, 
                           real colvector dnu, real colvector dlnlambda, real colvector dlnx, 
													 real colvector d2nu2, real colvector d2nulnlambda, real colvector d2nulnx, real colvector d2lnlambda2, real colvector d2lnlambdalnx, real colvector d2lnx2) {

	real colvector xl, power, diffuse, c, blockSize, mstar, mstarfrac, m, lnm, lnxlv, Slnxl, logRatios, cumLogRatios, lambdav, xv, lnKv, lnJv, lnpowersum, zeros, masspoints, tc, bigarg, bigord, big, dlnJlnlambda, d2lnJlnlambda2, d2lnJlnlambdalnx, d2lnJlnx2, dlnJlnx
	complex colvector termSums, _dnu, _dlnx
	complex matrix dTermSums, ddTermSums
	complex scalar lnBaseTermc, updatec, t2
	real scalar M, _blockSize, _c, _mstar, anybigarg, anybigord, x, lambda, xlambda, lnx, lnlambda, lnnu, lnxl, t, i, N,  lnBaseTerm, lnPeakTerm, LaplaceEpsilon, update, lnepsilon, todo, p, dpda, dpdx, _dnui, _dlnxi, d2lnKnulnlambda, d2lnKnulnx, nu, mnu, lnJ, lnK, tm, r, neghalfinvsqrtxl, dtdnu, dt2dnu2, twomm1, fournu2
	pointer (real colvector) scalar pSlnxl, pnuv, pdlnKnu, pdlnKlnlambda, pdlnKlnx

	N = max((rows(lnxv), rows(lnlambdav), rows(_nuv)))

	todo   =   _todo==.? 0 : _todo
	pnuv = reflect? &_nuv : &(-_nuv)

	lnepsilon = ln(epsilon(1)); LaplaceEpsilon = abs(invnormal(epsilon(.1)))

	lambdav = asdfExpNegInfty(lnlambdav)
	xv      = asdfExpNegInfty(lnxv)
	lnxlv = lnxv :+ lnlambdav; xl = xv :* lambdav

	c = ceil(1.5 :- *pnuv)  // index of first power series term with positive digamma(m+nu+1)^2-trigamma(m+nu+1) (and digamma(m+nu+1) and gamma(m+nu+1))
	mstar = -.5 * *pnuv :+ (sqrt(0.25 :* *pnuv :* *pnuv :+ xl))  // index of max power series term, anchor for expansion in both directions.
	t = floor(mstar); mstarfrac = mstar - t; mstar = t  // track fractional part of otherwise discrete mstar variable, for smoothing
	blockSize = ceil(LaplaceEpsilon * sqrt(mstar * .5))  // discretized too--but no smoothing for that yet!

	bigarg = blockSize :> 1e6  // use large-argument asymptotic or ascending power series expansion?
	power = 1 :- (big = bigarg)  // use asymptotic or ascending power series expansion?
	anybigord = any(bigord = bigarg :& *pnuv :/ xl :> 1e6)  // use large-order approximation?
	anybigarg = any(bigarg = bigarg :& !bigord)  // use large-argument asymptotic or ascending power series expansion?
	diffuse = power :& mstar :> 1e5  // many small power series terms vs few large ones: in latter case, for precision, sequence outer multipliers lnJ, lnK into computing individual terms before summing

	// factors outside sum; for precision break into two components with different scales: J should have scale of sum
	lnJv = xv :+ lambdav; if (anybigarg) lnJv = lnJv :- bigarg :* 2 * sqrt(xl)  // store lnJ negated to drop an operation from this line
	lnKv = anybigarg | anybigord? big   :* *pnuv*.5:*(lnxv:-lnlambdav) :- bigarg :* (.25*lnxlv :+ 1.43f89a3f0edd6X+000 /*ln(4*pi())*.5*/) +
	                              power :*  *pnuv :* (reflect? lnxv : lnlambdav) :
								                          *pnuv :* (reflect? lnxv : lnlambdav)
	if (todo) {
		dnu = dlnx = J(N, 1, 0)
		if (todo > 1) {
			d2nu2 = d2lnx2 = dnu
			d2nulnx = anybigord? .5 * bigord : dnu
		}

		/*dlnJdnu = 0*/
		dlnJlnlambda = anybigarg? (t=bigarg :* sqrt(lambdav:*xv)) :- lambdav : -lambdav
		dlnJlnx      = anybigarg?  t                              :-      xv : -     xv
		if (todo > 1) {
			// d2lnJnu2 = d2lnJnulnlambda = d2lnJnulnx = 0
			d2lnJlnlambda2   = anybigarg? (t = t*.5) :- lambdav : -lambdav
			d2lnJlnlambdalnx = anybigarg? (t       )            : 0
			d2lnJlnx2        = anybigarg? (t       ) :-      xv : -     xv
		}

		if (anybigarg | anybigord) {
			pdlnKnu       = &(bigarg :* (lnxv :- lnlambdav)*.5 :+ power :* (reflect? lnxv : lnlambdav))
			pdlnKlnlambda = &(big :* (t = -*pnuv*.5 :- .25*bigarg)); if (reflect==0) pdlnKlnlambda = &(*pdlnKlnlambda :+ power :* *pnuv)
			pdlnKlnx      = &(big :*  t :+ *pnuv                  ); if (reflect   ) pdlnKlnx      = &(*pdlnKlnx      :+ power :* *pnuv)
			if (todo > 1) {
				// d2lnKnu2 = 0
				d2lnKnulnlambda = power * (1 - reflect) - (t=big * .5)
				d2lnKnulnx      = power *      reflect  +  t
				// d2lnKlnlambda2 = 0
				// d2lnKlnlambdalnx = 0
				// d2lnKlnx2 = 0
			}
		} else {
			pdlnKnu       = reflect? &lnxv : &lnlambdav
			pdlnKlnlambda = reflect? &0    : pnuv
			pdlnKlnx      = reflect? pnuv  : &0
			if (todo > 1) {
				// d2lnKnu2 = 0
				d2lnKnulnlambda = 1 - reflect
				d2lnKnulnx      =     reflect
				// d2lnKlnlambda2 = 0
				// d2lnKlnlambdalnx = 0
				// d2lnKlnx2 = 0
			}
		}
	}

	if (lnxv != .) {
		lnx = lnxv[1]
    lnlambda = lnlambdav[1]
    lnxl= lnxlv[1]
    x = xv[1]
    lambda = lambdav[1]
    nu = (*pnuv)[1]
    _c = max((0, c[1]))
    lnJ = lnJv[1]
    lnK = lnKv[1]
		lnpowersum = J(N, 1, 0)

		for (i=N;i;i--) {  // main loop over observations
			if (rows(lnxv) > 1) {
				lnx = lnxv[i]
				x   = xv[i]
			}
			if (rows(lnlambdav) > 1) {
				lnlambda = lnlambdav[i]
				lambda   = lambdav  [i]
			}
			if (rows(lnxlv) > 1) {
				lnxl = lnxlv[i]
				lnJ  = lnJv [i]
			}
			if (rows(lnKv) > 1)
				lnK = lnKv[i]
			if (rows(*pnuv) > 1) {
				nu = (*pnuv)[i]
				_c = max((0, c[i]))
			}

			if (lnx == .) continue  // (deal with x = 0/lnx = . after loop)

			xlambda = x * lambda
			if (xlambda==. | nu==. | (reflect==0 & nu<0)) {
				lnpowersum[i] = .
				continue
			}

			if (lnlambda == . & (reflect | nu==0 | nu==-1)) {  // special case: central chi2, absorptive=reflective model for nu=0,-1
				lnpowersum[i] = -lngamma(nu + 1)  // sole survivor from log sum, in lambda->0 limit

				if (todo) {
					dnu[i] = -digamma(nu+1)  // and by default dlnx[i] = 0--because of lnlambda parameterization, can't really capture derivative wrt x, lambda when lnlambda=-infty
					if (todo > 1)
						d2nu2[i] = -trigamma(nu+1)
				}
				continue
			}

			if (bigord[i]) {  // large-order approximation, https://dlmf.nist.gov/10.41#i
				lnnu = ln(nu)
				lnpowersum[i] = nu + nu*(t = .5*lnxl - lnnu) - .5* lnnu - 1.d67f1c864beb4X-001 /*ln(sqrt(tau))*/
				if (todo) {
					dnu [i] = t - .5 / nu
					dlnx[i] = nu * .5
					if (todo > 1)
						d2nu2[i] = .5 / (nu*nu) - 1/nu
						// d2nulnx[i] = .5
						// d2lnx2 [i] = 0
				}
				continue
			}

			if (bigarg[i]) {  // if we'd need to sum more than ~2 million terms, use large-argument asymptotic expansion dlmf.nist.gov/10.40#E1
				termSums = lnBaseTermc = lnK - lnJ  // first summand is 1 (0 in logs)
        fournu2 = nu + nu; fournu2 = fournu2 * fournu2
				if (todo) {
        	twomm1 = tm = 1  // 2*m-1; first summand
          dtdnu = 0
          neghalfinvsqrtxl = -.5 / sqrt(xlambda)
					dTermSums = J(0, 2, C(.))  // since first summand (leaving aside J, K) is constant, logs of its derivatives are -infty
					if (todo > 1) {
          	dt2dnu2 = 0
						ddTermSums = J(0, 3, C(.))
          }
				}
				Slnxl = lnxl*.5 + 1.0a2b23f3bab73X+001 /*ln(8)*/
				m = 1

				do {
          updatec = ln(C((t = twomm1*twomm1 - fournu2)/(m+m))) - Slnxl
          if (Re(updatec) > 0) {  // asymptotic series diverged before converging--no good; give up!
            lnKv[i] = .  // force return of missing value
            break
          }
 					termSums = (lnBaseTermc = lnBaseTermc + updatec) \
					           termSums

					if (todo) {
          	r = t / (-8 * m)
						if (todo > 1)
              dt2dnu2 = neghalfinvsqrtxl * (r * dt2dnu2 + ((nu+nu)*dtdnu + tm) / m)

            dtdnu = neghalfinvsqrtxl * (r * dtdnu + nu/m * tm)
            tm = tm * r * neghalfinvsqrtxl
            
						dTermSums = ln(C(dtdnu)), lnBaseTermc + ln(C(m*-.5)) \ 
											  dTermSums

            if (todo > 1) {
            	t2 = ln(C(m*-.5))
              ddTermSums = (ln(C(dt2dnu2)), ln(C(tm))+t2, lnBaseTermc+t2+t2) \
													 ddTermSums
            }
					}
					if (lnBaseTermc==. | Re(lnBaseTermc) < lnepsilon)  // break if last update was small
						break
					++m; twomm1 = twomm1 + 2
				} while (1)

			} else {  // otherwise use ascending power series

        if (_c) {  // first few terms may pass negative arguments to gamma(m+nu+1); handle separately with code that is slower to allow for complex numbers in logs
					lnBaseTermc = (_c-1) * lnxl - lngamma(_c + nu) - lngamma(_c); if (diffuse[i]==0) lnBaseTermc = lnBaseTermc - lnJ + lnK  // term c-1; control order of operations in hope of improving precision
					termSums = lnBaseTermc
					if (todo) {
						_dnu  = lnBaseTermc + ln(C(t=-digamma(_c + nu)))  // logs of derivatives of base term, not log base term
						_dlnx = lnBaseTermc + (lnm=ln(_c-1))
						dTermSums = _dnu, _dlnx
						if (todo > 1)
							ddTermSums = lnBaseTermc + ln(C(t*t - trigamma(_c+nu))),
														_dnu       + lnm,
														_dlnx      + lnm
					}

					if (_c > 1) {  // compute additional negative terms if more than 1. Same logic as below, but slower because done in complex instead of real, to accomodate logs of negative ratios
						m = _c-1 :: 1; mnu = _c-1+nu :: 1+nu-.5
						logRatios = (lnm=ln(m)) + ln(mnu)
						cumLogRatios = quadrunningsum(logRatios, 1) - (1::_c-1) * lnxl
						termSums = lnBaseTermc + asdfLogSumExp(cumLogRatios) \ termSums
						if (todo) {
							_dnu  = cumLogRatios +      ln(C(t=digamma(mnu)))
							_dlnx = cumLogRatios + (lnm = _c==2? . : lnm[|2\.|] \ .)
							dTermSums = lnBaseTermc :+ (1.921fb54442d18X+001i /*pi i*/ + asdfLogSumExp(_dnu )  ,  // logs of nu derivatives
													                                                 asdfLogSumExp(_dlnx)) \  // logs of x/lambda derivatives
													dTermSums
							if (todo > 1)
								ddTermSums = lnBaseTermc :+ (                                 asdfLogSumExp(cumLogRatios + ln(C(t:*t - trigamma(mnu))))  ,
								                             1.921fb54442d18X+001i /*pi i*/ + asdfLogSumExp(_dnu         + lnm                        )  ,
								                                                              asdfLogSumExp(_dlnx        + lnm                        )) \
                             ddTermSums
						}
					}
				} else {
					termSums = J(0, 1, 0)
					if (todo) {
						dTermSums = J(0, 2, C(0))
						if (todo > 1)
							ddTermSums = J(0, 3, 0)
					}
				}

				_mstar = max((_c, mstar[i]))
				lnPeakTerm = _mstar<1e5? _mstar *  lnxl - (lngamma(_mstar + nu + 1) + lngamma(_mstar + 1)) :  // use Stirling's approx for large m* because terms mostly cancel out
				                         _mstar * (lnxl - (ln     (_mstar+nu      ) + ln     (_mstar    ))) - 1.d67f1c864beb4X+000 /*ln(2*pi())*/  + nu -.5*ln(_mstar) + ( -(nu+.5)*ln(_mstar+nu)) + (_mstar+_mstar)
				if (diffuse[i]==0) lnPeakTerm = lnPeakTerm - lnJ + lnK
				termSums = lnPeakTerm \ termSums

				if (todo) {
					_dnu  = lnPeakTerm + ln(C(t=-digamma(_mstar + nu + 1)))  // logs of derivatives of peak term, not log peak term
					_dlnx = lnPeakTerm + (lnm=ln(_mstar))
					dTermSums = dTermSums \ (_dnu, _dlnx)
					if (todo > 1)
						ddTermSums = ddTermSums \ (lnPeakTerm + ln(C(t*t - trigamma(_mstar+nu+1))),
												               _dnu       + lnm,
												               _dlnx      + lnm)
				}

				_blockSize = max((100, blockSize[i]))
				Slnxl = (1::_blockSize) * lnxl
				if (_mstar > _c) {
					m = _mstar :: 1 + max((_c, _mstar-_blockSize))  // block of indexes for positive terms preceding peak

					lnBaseTerm = lnPeakTerm
					do {
            M = rows(m)
            if (m[M]<epsilon(nu))
              mnu = J(M,1,nu)
            else {
              mnu = m[1]+nu :: m[M]+nu-.5  // faster than m:+nu
              if (length(mnu) > M) mnu = mnu[|.\M|]
            }
						pSlnxl = &(M==_blockSize? Slnxl : Slnxl[|.\M|])
						logRatios = (lnm=ln(m)) + ln(mnu)  // log of ratio of successive denominators in power series
						cumLogRatios = quadrunningsum(logRatios) - *pSlnxl
						termSums = (update = lnBaseTerm + asdfLogSumExp(cumLogRatios)) \ termSums

						if (todo) {
							_dnu  = cumLogRatios + ln(t=digamma(mnu))
							_dlnx = cumLogRatios + (lnm = M==1? ln(m-1) : lnm[|2\.|]\ln(m[M]-1))
							dTermSums = lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/ + asdfLogSumExp(_dnu ) ,  // logs of nu derivatives of sum of terms
													lnBaseTerm                                  + asdfLogSumExp(_dlnx) \  // logs of x/lambda derivatives
                          dTermSums
							if (todo > 1)
								ddTermSums = lnBaseTerm                                 + asdfLogSumExp(cumLogRatios + ln(t:*t - trigamma(mnu))) ,
														 lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/+ asdfLogSumExp(_dnu         + lnm                     ) ,
														 lnBaseTerm                                 + asdfLogSumExp(_dlnx        + lnm                     ) \
                             ddTermSums
						}
            
            if (m[M]==_c+1 | update==. | rows(m)==1) break  // reached early, negative terms, or hit stopping condition and did last fractional term for smoothing

						m = lnepsilon > cumLogRatios[1]+(lnBaseTerm - colmax(Re(termSums))) + ln(_blockSize)? m[M]-1 : m[M]-1::max((m[M]-_blockSize , _c+1)) // if hit stopping condition, do one last, fractionally weighted term
            
            lnBaseTerm = m[1]<1e5? m[1] *  lnxl - lngamma(m[1] + nu + 1) - lngamma(m[1] + 1) :  // use Stirling's approx for large m because terms mostly cancel out
						                       m[1] * (lnxl -(ln     (m[1] + nu    ) + ln     (m[1]    ))) - 1.d67f1c864beb4X+000 /*ln(2*pi())*/  +nu -.5*ln(m[1]) + ( -(nu+.5)*ln(m[1]+nu)) + (m[1]+ m[1])
						if (diffuse[i]==0) lnBaseTerm = lnBaseTerm - lnJ + lnK
            if (rows(m)==1) lnBaseTerm = lnBaseTerm + ln(1-mstarfrac[i])  // fractionally weighted, last singleton term block
					} while (1)
				}

				m = _mstar+1 :: _mstar+_blockSize  // block of indexes for terms following peak
				lnBaseTerm = lnPeakTerm
				do {
          M = rows(m)
					if (m[1]<epsilon(nu))
            mnu =  J(M,1,nu)
          else {
          	mnu = m[1]+nu :: m[M]+nu+.5  // faster than m:+nu
            if (length(mnu) > M) mnu = mnu[|.\M|]
          }
					pSlnxl = &(M==_blockSize? Slnxl : Slnxl[|.\M|])
					logRatios = (lnm=ln(m)) + ln(mnu)  // log of ratio of successive denominators in power series
					cumLogRatios = *pSlnxl - quadrunningsum(logRatios)
					termSums = (update = lnBaseTerm + asdfLogSumExp(cumLogRatios)) \ termSums
					if (todo) {
						if (m[M]>=epsilon(nu)) {
              mnu = m[1]+nu+1 :: m[M]+nu+1.5
              if (length(mnu) > M) mnu = mnu[|.\M|]
            }
						_dnu = cumLogRatios + ln(t=digamma(mnu))
						_dlnx   = cumLogRatios + lnm
						dTermSums = lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/ + asdfLogSumExp(_dnu ) ,  // logs of nu derivatives
												lnBaseTerm                                  + asdfLogSumExp(_dlnx) \  // logs of x/lambda derivatives
                        dTermSums
							if (todo > 1)
								ddTermSums = lnBaseTerm +                                  asdfLogSumExp(cumLogRatios + ln(t:*t - trigamma(mnu))) ,
														 lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/ + asdfLogSumExp(_dnu         + lnm                      ) ,
														 lnBaseTerm +                                  asdfLogSumExp(_dlnx        + lnm                      ) \
                             ddTermSums
					}

          if (update==. | rows(m)==1) break  // hit stopping condition and did last fractional term for smoothing

					m = lnepsilon > cumLogRatios[_blockSize]+(lnBaseTerm-colmax(Re(termSums)))+ln(_blockSize)? m[1]+_blockSize : m[1]+_blockSize::m[M]+_blockSize  // if hit stopping condition, do one last, fractionally weighted term

					lnBaseTerm = m[1]<1e5? (m[1]-1) *  lnxl - lngamma(m[1]     + nu) - lngamma(m[1]  ) :  // use Stirling's approx for large m* because terms mostly cancel out
					                       (m[1]-1) * (lnxl -(ln     (m[1] - 1 + nu) + ln     (m[1]-1))) - 1.d67f1c864beb4X+000 /*ln(2*pi())*/  +nu -.5*ln(m[1]-1) + ( -(nu+.5)*ln(m[1]-1+nu)) + (m[1]+m[1]-2)
					if (diffuse[i]==0) lnBaseTerm = lnBaseTerm - lnJ + lnK
          if (rows(m)==1) lnBaseTerm = lnBaseTerm + ln(mstarfrac[i])  // fractionally weighted, last singleton term block
				} while (1)
			}  // end ascending power series

			lnpowersum[i] = eltype(termSums)=="complex"? Re(asdfLogSumExp(termSums)) : asdfLogSumExp(termSums)
      t = bigarg[i]? 0 : lnpowersum[i]
			if (todo) {
				tc = Re(exp(asdfLogSumExp(dTermSums) :- t))  // log(deriv of log sum) = log(deriv of sum) - log sum
				dnu [i] = _dnui  = tc[1]
				dlnx[i] = _dlnxi = tc[2]
				if (todo > 1) {
					tc = Re(exp(asdfLogSumExp(ddTermSums) :- t))
					d2nu2  [i] = tc[1] - _dnui  * _dnui
					d2nulnx[i] = tc[2] - _dnui  * _dlnxi
					d2lnx2 [i] = tc[3] - _dlnxi * _dlnxi
				}
			}
		}

		if (todo) {
			dnu       = *pdlnKnu       :+/*dlnJdnu=0    :+*/ dnu  ; if (reflect==0) dnu = -dnu  // adjust for initial negating of nu 
			dlnlambda = *pdlnKlnlambda :+ (dlnJlnlambda :+   dlnx)  // add dlnJ first for precision (hoping)
			dlnx      = *pdlnKlnx      :+ (dlnJlnx      :+   dlnx)  // derivatives of sum are same w.r.t lnlambda and lnx
			if (todo > 1) {
			//d2nu2         =                                        d2nu2         
				d2nulnlambda  =  d2lnKnulnlambda                  :+   d2nulnx 
				d2nulnx       =  d2lnKnulnx                       :+   d2nulnx
				d2lnlambda2   =                  d2lnJlnlambda2   :+   d2lnx2
				d2lnlambdalnx =                  d2lnJlnlambdalnx :+   d2lnx2
				d2lnx2        =                  d2lnJlnx2        :+   d2lnx2
				if (reflect==0) {
					d2nulnlambda = -d2nulnlambda; d2nulnx = -d2nulnx
				}
			}
		}
	}

	lnJv = lnJv :* diffuse; lnKv = lnKv :* diffuse

	zeros = lnxv:==.
	masspoints = zeros :& _nuv :< 0 :& ((reflect==0) :| round(_nuv):==(_nuv))
	zeros = zeros :& !masspoints
	if (zeros==1) {  // if len(zeros)=len(masspoints)=max(len(nu),len(X))==1 then *all* (lambda-varying) obs have x=0
		zeros = .; masspoints = J(0,1,0)  // as a matrix index, . means all rows
	} else {
		zeros = selectindex(zeros)
		if (cols(zeros)==0) zeros = J(0,1,0)
		masspoints = selectindex(masspoints)
		if (cols(masspoints)==0) masspoints = J(0,1,0)
	}
	N = zeros==.? rows(lnKv) : rows(zeros)

	if (N) {  // x=0, but no mass point there
		t = J(N, 1, .)
		lnKv[zeros] = t
		if (todo) {
			dnu[zeros] = dlnlambda[zeros] = t
			if (todo > 1)
				d2nulnx[zeros] = d2lnlambdalnx[zeros] = d2lnx2[zeros] = d2nu2[zeros] = d2nulnlambda[zeros] = d2lnlambda2[zeros] = t
		}
	}

	N = masspoints==.? rows(lnKv) : rows(masspoints)
	if (N) {
		if (rows(lambdav) == 1) {
			p = gammaptail(-_nuv, lambdav)
			lnKv[masspoints] = J(N, 1, ln(p))  //p<.5? ln(p) : ln1m(gammap(-_nuv, lambdav))  would maximize precision...
			lnJv[masspoints] = J(N, 1, 0)
		} else {
			lambdav = lambdav[masspoints]
			p = gammaptail(-_nuv, lambdav)
			lnKv[masspoints] = ln(p)  //p<.5? ln(p) : ln1m(gammap(-_nuv, lambdav))  would maximize precision...
      lnJv[masspoints] = J(N, 1, 0)
		}
		if (todo) {
			t = J(N, 1, .)
			dnu      [masspoints] = _dnui  =  reflect? t :  (dpda = dgammapda(-_nuv, lambdav)) :/ p    // in reflective (noncentral chi2) case, nu=-1 and p is discontinuous w.r.t. nu at x=0
			dlnlambda[masspoints] = _dlnxi = -lambdav    :* (dpdx = dgammapdx(-_nuv, lambdav)) :/ p
			if (todo > 1) {
				d2nulnx     [masspoints] = d2lnlambdalnx[masspoints] = d2lnx2[masspoints] = t
				d2nu2       [masspoints] = reflect? t : _dnui  :* (               -dgammapdada(-_nuv, lambdav):/dpda - _dnui )
				d2nulnlambda[masspoints] = reflect? t : _dnui  :* (     lambdav :* dgammapdadx(-_nuv, lambdav):/dpda - _dlnxi)
				d2lnlambda2 [masspoints] =              _dlnxi :* (1 :+ lambdav :* dgammapdxdx(-_nuv, lambdav):/dpdx - _dlnxi)
			}
		}
	}
	return(rows(lnpowersum)? lnKv + (lnpowersum :- lnJv) : lnKv :- lnJv)
}

// causes large-argument series to instantly diverge, but seems to require impractically many terms for power series approach, and yet order is small relative to  argument, so seemingly large-order approx is bad too:
//lnPDFFeller( 4621016466 ,  27.26823829  , 27.27523152)

// mimic Stata's nchi2den()
real colvector mynchi2den(real scalar k, real colvector lambda, real colvector x)  // should match Stata's nchi2den()
	return (exp(-1.62e42fefa39efX-001 :+ lnPDFFeller(k*.5-1, ln(lambda):-1.62e42fefa39efX-001, ln(x):-1.62e42fefa39efX-001 /*ln 1/2*/)))


// crude power series implementation of nchi2 CDF/Marcum Q function, with extension to Feller distribution
real colvector CDFFeller(real colvector nuv, real colvector lnlambdav, real colvector lnxv, | real scalar _reflect) {
	real scalar i, N, x, reflect, nu, lambda, lnlambda; real colvector m, lambdav, xv, retval
	
	N = max((rows(nuv), rows(lnxv), rows(lnlambdav)))
	reflect = editmissing(_reflect, 1)

	m = 1e6::0

	lambdav = asdfExpNegInfty(lnlambdav); xv = asdfExpNegInfty(lnxv)
	lambda = lambdav[1]; lnlambda = lnlambdav[1]; x = xv[1]; nu = nuv[1]
	retval = J(N, 1, .)
	for (i=N;i;i--) {
		if (rows(lnxv) > 1)
			x = xv[i]
		if (rows(lnlambdav) > 1) {
			lambda = lambdav[i]
			lnlambda = lnlambdav[i]
		}
		if (rows(nuv)>1)
			nu = nuv[i]
		if (reflect & nu<=-1)
			continue

		retval[i] = quadcolsum(exp((reflect? m      *lnlambda - lngamma(m :+  1    ) + ln(gammap((nu+1) :+ m, x))  :
		                                     (m:-nu)*lnlambda - lngamma(m :+ (1-nu)) + ln(gammap(    1  :+ m, x))) :- lambda))

		if (nu < 0 & (reflect==0 | round(nu)==(nu)))  // mass at 0?
			retval[i] = retval[i] + gammaptail(-nu, lambda)
	}
	return(retval)
}

// returns wdith-2 matrix, col 1 for mean, col 2 for variance, of Feller/noncentral chi2 distribution
real matrix meanVarianceFeller(real colvector nuv, real colvector lnlambdav, | real scalar reflect) {
	real colvector lambdav, F, f, mean, lnu
	lambdav = exp(lnlambdav)
	lnu = lambdav :+ nuv
	if (reflect == 0) {
		f = lambdav :* gammaden(-nuv, 1, 0, lambdav)
		F =            gammap  (-nuv,       lambdav)
		mean = f + (lnu :+ 1) :* F
		return ((mean , (lnu :+ 3) :* f + (lambdav :+ (nuv :+ 1) :* (lnu :+ 2)) :* F - mean:*mean))
	}
	mean = lnu:+1
	return ((mean , mean :+ lambdav))
}

// use secant search to invert Feller/noncentral chi2 distribution CDF function at percentiles p
real matrix invCDFFeller(real vector p, real colvector nuv, real colvector lnlambdav, real scalar reflect) {
	real scalar _lo, _hi, i, j, N, M, masspoint, mid, p_mid, _p, p_lo, p_hi, nu, lnlambda; real matrix q; real colvector mv, z, lo, hi

	N = max((rows(nuv), rows(lnlambdav))); M = length(p)
	q = J(N, M, .)
	mv = meanVarianceFeller(nuv, lnlambdav, reflect)
	z = invnormal(rows(p)==1? p' : p) * mv[,2]
	nu = nuv[1]; lnlambda = lnlambdav[1]
	for (i=N; i; i--) {
		if (rows(nuv      ) > 1) nu       = nuv      [i]
		if (rows(lnlambdav) > 1) lnlambda = lnlambdav[i]
		masspoint = CDFFeller(nu, lnlambda, ., reflect)
		lo = editmissing(ln(mv[i,1] :- z), ln(1e-6))
		hi = editmissing(ln(mv[i,1] :+ z), ln(1e-6))
		for (j=M;j;j--)
			if ((_p = p[j]) > masspoint) {  // default to 0 if in mass point
				_lo = lo[j]; _hi = hi[j]
				while (_p < (p_lo = CDFFeller(nu, lnlambda, _lo, reflect)))
					_lo = _lo - ln(10)
				while (_p > (p_hi = CDFFeller(nu, lnlambda, _hi, reflect)))
					_hi = _hi + ln(10)

					if (_lo==. | _hi==.)
            q[i,j] = .  // can't bracket quantile with double-precision numbers
				else {
					do {
						mid = _lo + (_p - p_lo) / (p_hi - p_lo) * (_hi - _lo)  // interpolate linearly
						if (mreldif(_lo, mid)<1e-12 | mreldif(_hi, mid)<1e-12)
							break

						p_mid = CDFFeller(nu, lnlambda, mid, reflect)
						if ((p_mid < _p) == (p_lo < _p)) {
							 _lo =   mid
							p_lo = p_mid
						} else {
							 _hi =   mid
							p_hi = p_mid
						}
					} while (1)
					q[i,j] = mid
				}
			}
	}
	return(q)
}


// fit Y = ((Y0^-B + s/delta)*exp(-delta*Bt) - s/delta)^-(1/B), minimizing sum squared errors in logs; Y0 is a parameter too
void BernoulliNLSStatic(real scalar todo, real rowvector params, real colvector lny, real colvector t, real colvector wt, real scalar sig2, real matrix g, real matrix H) {
	real colvector e, ewt, p, edivp, lnp, mdeltaBt, dY0, expmdeltaBt; real scalar s, B, delta, Y0, sdivdelta, ds, Y0B
	pragma unset H
	s = params[1]; B = params[2]; delta = params[3]; Y0 = params[4]
	sdivdelta = s / delta
	mdeltaBt = (-delta) * B * t
	expmdeltaBt = exp(mdeltaBt)
	Y0B = Y0^-B
	p = (Y0B + sdivdelta) * expmdeltaBt :- sdivdelta
	lnp = ln(p)
	e = lny + lnp/B
	ewt = e :* wt
	sig2 = e ' ewt
	if (todo) {
		edivp = ewt :/ p
		ds = edivp ' exp(mdeltaBt-1) / (delta*B)
		dY0 = edivp :* expmdeltaBt
		g = ds, 
		    -(dY0 ' (Y0B*ln(Y0):+(Y0B*delta+s)*t) + ewt ' lnp / B) / B, 
				-sdivdelta*ds - (Y0B + sdivdelta) * (dY0 ' t), 
				-Y0B/Y0 * quadcolsum(dY0, 1)
		g = g + g																																																																																																																		
	}
}

mata mlib create lasdf, dir("`c(sysdir_plus)'l") replace
mata mlib add lasdf *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
