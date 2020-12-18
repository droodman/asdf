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
void asdfEst::setData(real colvector _Y0, real colvector _Y, real colvector _tDelta) {
	lnY0 = ln(_Y0)
	lnY  = ln(_Y )
	zerofill = J(length(zeroind = selectindex(_Y:==0)), 1, 0)  // not used in GBM
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
		if (length(ind = selectindex(*pb :== 0)))
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

real scalar asdfEststickyfeller::getlf() return(0)  // only an lf0 estimator

// _Y0, _Y, _tDelta = inital & final values, and time gaps
void asdfEststickyfeller::setData(real colvector Y0, real colvector Y, real colvector tDelta) {
	pY0 = &Y0; pY  = &Y; ptDelta = &tDelta
	S.setData(tDelta, Y0, Y)
}

string rowvector asdfEststickyfeller::getParamEstNames()
	return (("lna", "b", "nu", "lnmu"))

string rowvector asdfEststickyfeller::getParamFixedNames()
	return ("")

void asdfEststickyfeller::setFixedParam(string scalar name, numeric value) {}

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
	lnf = S.lnStickyFeller(exp(*plna), *plna, *pb, *pnu, exp(*plnmu), *plnmu)
}


//
// class for sticky squared Bessel, as special case of sticky Feller
//

class asdfEststickysqbessel extends asdfEststickyfeller {
	void new()
	virtual void getScalarParams(), getSmatrixParams()
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


// lf2 likelihood evaluator for all estimation classes
// 1st argument is a moptimize() obect whose 1st userinfo item is the asdfEst obect
function asdflf2(transmorphic scalar M, real scalar todo, real rowvector b, real colvector lnf, real matrix g, real matrix H) {
	class asdfEst scalar S; struct smatrix matrix h, params; real scalar d, i, j; real matrix t
	pragma unset h
	S = moptimize_util_userinfo(M, 1)
	d = rows(params = S.getParams())
	for (i=d; i; i--)
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

// sum columns of potentially negative numbers stored in logs, avoiding overflow, treating missing as log of 0
numeric rowvector asdfquadLogSumExp(numeric matrix x) {
	real rowvector shift
	if (rows(x)==0) return(J(1,cols(x),0))
	if (rows(x)==1) return(x)
	shift = ln(maxdouble()/rows(x)) :- colmax(eltype(x)=="real"? x : Re(x))
//	shift = shift - (shift:>0):*shift  // only downshift, to present overflow; shifting can prevent underflow & overflow but can also reduce precision if the shifter is much larger than entries
	return (ln(quadcolsum(exp(x :+ shift))) - shift)
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
		lnx = lnxv[1]; lnlambda = lnlambdav[1]; lnxl= lnxlv[1]; x = xv[1]; lambda = lambdav[1]; nu = (*pnuv)[1]; _c = max((0, c[1])); lnJ = lnJv[1]; lnK = lnKv[1]
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
			if (xlambda==. | nu==. | (!reflect & nu<0)) {
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
						termSums = lnBaseTermc + asdfquadLogSumExp(cumLogRatios) \ termSums
						if (todo) {
							_dnu  = cumLogRatios +      ln(C(t=digamma(mnu)))
							_dlnx = cumLogRatios + (lnm = _c==2? . : lnm[|2\.|] \ .)
							dTermSums = lnBaseTermc :+ (1.921fb54442d18X+001i /*pi i*/ + asdfquadLogSumExp(_dnu )  ,  // logs of nu derivatives
													                                                 asdfquadLogSumExp(_dlnx)) \  // logs of x/lambda derivatives
													dTermSums
							if (todo > 1)
								ddTermSums = lnBaseTermc :+ (                                 asdfquadLogSumExp(cumLogRatios + ln(C(t:*t - trigamma(mnu))))  ,
								                             1.921fb54442d18X+001i /*pi i*/ + asdfquadLogSumExp(_dnu         + lnm                        )  ,
								                                                              asdfquadLogSumExp(_dlnx        + lnm                        )) \
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
						termSums = (update = lnBaseTerm + asdfquadLogSumExp(cumLogRatios)) \ termSums

						if (todo) {
							_dnu  = cumLogRatios + ln(t=digamma(mnu))
							_dlnx = cumLogRatios + (lnm = M==1? ln(m-1) : lnm[|2\.|]\ln(m[M]-1))
							dTermSums = lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/ + asdfquadLogSumExp(_dnu ) ,  // logs of nu derivatives of sum of terms
													lnBaseTerm                                  + asdfquadLogSumExp(_dlnx) \  // logs of x/lambda derivatives
                          dTermSums
							if (todo > 1)
								ddTermSums = lnBaseTerm                                 + asdfquadLogSumExp(cumLogRatios + ln(t:*t - trigamma(mnu))) ,
														 lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/+ asdfquadLogSumExp(_dnu         + lnm                     ) ,
														 lnBaseTerm                                 + asdfquadLogSumExp(_dlnx        + lnm                     ) \
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
					termSums = (update = lnBaseTerm + asdfquadLogSumExp(cumLogRatios)) \ termSums
					if (todo) {
						if (m[M]>=epsilon(nu)) {
              mnu = m[1]+nu+1 :: m[M]+nu+1.5
              if (length(mnu) > M) mnu = mnu[|.\M|]
            }
						_dnu = cumLogRatios + ln(t=digamma(mnu))
						_dlnx   = cumLogRatios + lnm
						dTermSums = lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/ + asdfquadLogSumExp(_dnu ) ,  // logs of nu derivatives
												lnBaseTerm                                  + asdfquadLogSumExp(_dlnx) \  // logs of x/lambda derivatives
                        dTermSums
							if (todo > 1)
								ddTermSums = lnBaseTerm +                                  asdfquadLogSumExp(cumLogRatios + ln(t:*t - trigamma(mnu))) ,
														 lnBaseTerm + 1.921fb54442d18X+001i /*pi i*/ + asdfquadLogSumExp(_dnu         + lnm                      ) ,
														 lnBaseTerm +                                  asdfquadLogSumExp(_dlnx        + lnm                      ) \
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

			lnpowersum[i] = eltype(termSums)=="complex"? Re(asdfquadLogSumExp(termSums)) : asdfquadLogSumExp(termSums)
      t = bigarg[i]? 0 : lnpowersum[i]
			if (todo) {
				tc = Re(exp(asdfquadLogSumExp(dTermSums) :- t))  // log(deriv of log sum) = log(deriv of sum) - log sum
				dnu [i] = _dnui  = tc[1]
				dlnx[i] = _dlnxi = tc[2]
				if (todo > 1) {
					tc = Re(exp(asdfquadLogSumExp(ddTermSums) :- t))
					d2nu2  [i] = tc[1] - _dnui  * _dnui
					d2nulnx[i] = tc[2] - _dnui  * _dlnxi
					d2lnx2 [i] = tc[3] - _dlnxi * _dlnxi
				}
			}
		}

		if (todo) {
			dnu       = *pdlnKnu       :+/*dlnJdnu=0    :+*/ dnu  ; if (!reflect) dnu = -dnu  // adjust for initial negating of nu 
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
	masspoints = zeros :& _nuv :< 0 :& (!reflect :| round(_nuv):==(_nuv))
	zeros = zeros :& masspoints:==0
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
lnPDFFeller( 4621016466 ,  27.26823829  , 27.27523152)

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

		if (nu < 0 & (!reflect | round(nu)==(nu)))  // mass at 0?
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
