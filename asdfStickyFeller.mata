global lnStickyFeller_M 60  // max term index in power series expansion of U()
global lnBesselKPower_M 22  // max term index in power series expansion of BesselK()
global lnUAsymParam_M 12 // max term index in large-parameter expansion of U(); 2S+3 = last Bernoulli number needed
global lnStickyFeller_N 28  // number of points sampled on parabola of integration

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

class clsStickyFeller {
	real colvector lnStickyFeller_m, lnStickyFeller_lnm, lnUAsymParam_B, lnUAsymParam_kB, lnBesselKPower_0toM, lnBesselKPower_2Mto0, lnBesselKAsym_k, lnBesselKAsym_twokm1sq
	complex rowvector lnStickyFeller_tlambdalnu, lnStickyFeller_u2
	pointer(real colvector) colvector lnUAsymParam_PT
	real rowvector lnUAsymParam_p, lnUAsymParam_q

	real colvector lnStickyFeller()
	complex matrix lnStickyGreen()
	complex scalar hypergeomdiff()
	complex rowvector lnUAsymParam(), lnBesselKAsym(), lnBesselKPower()
	void new()
}

void clsStickyFeller::new() {
	lnStickyFeller_m = 1::$lnStickyFeller_M; lnStickyFeller_lnm = ln(lnStickyFeller_m)

  complex rowvector u
	u = C(1, (.5..27.5) * 1.6c782f6914edbX-003 /*h*/)  // really (1 + ui)
	lnStickyFeller_u2 = u :* u
	lnStickyFeller_tlambdalnu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ * lnStickyFeller_u2 + ln(u) :- 1.16c0cdab51479X+001 /* + ln 2*2*h/2pi  tack this on here for efficiency */

  lnBesselKPower_0toM = 0::$lnBesselKPower_M  // stuff for power series
  lnBesselKPower_2Mto0 = lnBesselKPower_0toM * 2

	lnBesselKAsym_k = 0::7  // 7 comes from using 1st 8 terms of series
	lnBesselKAsym_twokm1sq = 2*lnBesselKAsym_k:-1; lnBesselKAsym_twokm1sq = lnBesselKAsym_twokm1sq :* lnBesselKAsym_twokm1sq
	lnBesselKAsym_k = ln(8 * lnBesselKAsym_k)  // 8 comes from asymptotic formula

	lnUAsymParam_p = lnUAsymParam_q = J($lnUAsymParam_M,1,1)

	// B_i / i!, i=1 to 99, B_i = Bernoulli numbers
	lnUAsymParam_B = -.5 \ 0.0833333333333333 \ 0 \ -.00138888888888889 \ 0 \ 3.30687830687831e-5 \ 0 \ -8.26719576719577e-7 \ 0 \ 2.08767569878681e-8 \ 0 \ -5.28419013868749e-10 \ 0 \ 1.33825365306847e-11 \ 0 \ -3.38968029632258e-13 \ 0 \ 8.58606205627785e-15 \ 0 \ -2.17486869855806e-16 \ 0 \ 5.50900282836023e-18 \ 0 \ -1.39544646858125e-19 \ 0 \ 3.53470703962947e-21 \ 0 \ -8.95351742703755e-23 \ 0 \ 2.26795245233768e-24 \ 0 \ -5.7447906688722e-26 \ 0 \ 1.45517247561486e-27 \ 0 \ -3.68599494066531e-29 \ 0 \ 9.33673425709505e-31 \ 0 \ -2.36502241570063e-32 \ 0 \ 5.99067176248213e-34 \ 0 \ -1.51745488446829e-35 \ 0 \ 3.84375812545419e-37 \ 0 \ -9.73635307264669e-39 \ 0 \ 2.46624704420068e-40 \ 0 \ -6.24707674182074e-42 \ 0 \ 1.58240302446449e-43 \ 0 \ -4.00827368594894e-45 \ 0 \ 1.01530758555696e-46 \ 0 \ -2.57180415824187e-48 \ 0 \ 6.51445603523381e-50 \ 0 \ -1.65013099068965e-51 \ 0 \ 4.17983062853948e-53 \ 0 \ -1.05876346677029e-54 \ 0 \ 2.68187919126077e-56 \ 0 \ -6.79327935110742e-58 \ 0 \ 1.72075776166814e-59 \ 0 \ -4.35873032934889e-61 \ 0 \ 1.10407929036847e-62 \ 0 \ -2.79666551337813e-64 \ 0 \ 7.08403650167947e-66 \ 0 \ -1.79440740828922e-67 \ 0 \ 4.5452870636111e-69 \ 0 \ -1.15133466319821e-70 \ 0 \ 2.91636477109236e-72 \ 0 \ -7.38723826349734e-74 \ 0 \ 1.8712093117638e-75 \ 0 \ -4.7398285577618e-77 \ 0 \ 1.20061259933545e-78 \ 0
	lnUAsymParam_kB = lnUAsymParam_B[|2\.|] :* (1::rows(lnUAsymParam_B)-1) \ 0

	// Pascal's triangle
	lnUAsymParam_PT = &(1\1)
	real matrix PT; real scalar s
	PT = comb(2..$lnUAsymParam_M, 0::$lnUAsymParam_M)
	for (s=1; s<=$lnUAsymParam_M-1; s++)
		lnUAsymParam_PT = lnUAsymParam_PT \ &(PT[|.,s\s+2,s|])
}


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007), Weideman (2010), parabolic contour method to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real colvector t, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu) {
  real colvector _mu, lnp0, negbt, lnattilde, lnm, lnx, lny, nulnx, nulny, zerofill, zeroind, missingfill, lnJacobian, bax, bay; real scalar ba

	if (nu < -1 | nu > 0 | mu < 0)
		return (.)

	zerofill = J(length(zeroind = selectindex(y:==0)), 1, 0)
	missingfill = J(length(zeroind), 1, .)
	lnx = ln(x); lny = ln(y)

	if (b)
		lnattilde = ln(a / (-b) * expm1(negbt = (-b) * t))
	else {
		lnattilde = ln(a * t)
		negbt = 0
	}
	lnJacobian =  negbt :- lnattilde

	if (mu == .) {
		lnp0 = lnPDFFeller(nu, lnx - lnattilde, lny + lnJacobian, 1)  // reflecting transition density
		if (length(zeroind)) lnJacobian[zeroind] = zerofill
		return (lnp0 + lnJacobian)
	}
	lnp0 = lnPDFFeller(nu, lnx - lnattilde, lny + lnJacobian, 0)  // absorbing transition density
	if (rows(zerofill)) lnJacobian[zeroind] = zerofill
	lnp0 = lnp0 + lnJacobian
	if (mu == 0 | nu == 0 | nu == -1)
		return (lnp0)
	if (rows(zerofill)) lnp0[zeroind] = missingfill  // if y=0, "absorbing term" drops out

	ba = abs(b / a); bax = ba *x; bay = ba * y
	nulnx = nu * lnx; nulny = nu * lny
	_editmissing(lnm = nulny - (b < 0? bay : bax) :- ln(a), -ln(mu))  // log speed measure; takes special value of 1/mu when y = 0; b>0 modification is a hack to multiply by exp(-b/a*(x+y))

/*	log_epsilon = ln(1e-15)  // Logic extracted from Garrappa's mlf routine for Matlab
	log_eps = ln(epsilon(1))
	w = sqrt(log_eps / (log_eps - log_epsilon))  // half-width of integration range, needed to assure given precision
	N = ceil(-w * log_epsilon / (2*pi())) + 1    // half the number of integration points
	h = w :/ (N-.5)                              // width of bars in Riemann integral; do an even number of integration points since only integrating over half of parabola
  u = C(1, (.5..27.5) * 1.6c782f6914edbX-003)  // really 1 + ui */

  _mu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ :/ t  // different "mu" -- as in Garrappa (2015), p. 1364.

	
	// Done computing asborbing term and prepping for parabolic path integral for sticky term. Now compute sticky term along parabola.

	real colvector tp1, tp2, denompx1, denompx2, denompy1, denompy2, numerp1, numerp2, s, xhi, yhi
	real scalar muj, _lnC, lnbax, lnbay, i, j, lngammanegnu, _nulnx, _nulny, beta, _bax, _bay, qa, qb, qc, qd, smaxx, smaxy, Ream1, ta, ta2, ta3, absa2, Ima2, _ix, _iy, absb, anypowerx, anyasymargx, anyasymparamx, anypowery, anyasymargy, anyasymparamy, powerx, asymargx, asymparamx, powery, asymargy, asymparamy
	complex matrix alpha
	complex scalar lnC, lnM, _alpha, _alpha1beta, lngammaalpha1beta, retval
	complex rowvector lnUAsymParamx, lnUAsymParamy
	complex colvector terms, numera

	beta = 1 + nu

	if ((absb=abs(b)) > 1e-10) {
		alpha = _mu * (lnStickyFeller_u2 / -b); if (b > 0) alpha = beta :- alpha  // faster to just compute alpha row by row? (but will still need vessel for retval)

		xhi = bax :> 14.18699
		yhi = bay :> 14.18699
		
		lngammanegnu = lngamma(-nu)
		_lnC = lngamma(beta) - ln(-nu) - lngammanegnu - nu * ln(abs(b)/a)

		tp1 = ln(beta   :: $lnStickyFeller_M.5 + nu) + lnStickyFeller_lnm  // faster than ln(nu :+ m)
		tp2 = ln(1 - nu :: $lnStickyFeller_M.5 - nu) + lnStickyFeller_lnm  // needed for power series only

		for (j=rows(x);j;j--) {
			lnbax = ln(_bax = bax[j]); _nulnx = nulnx[j]
			lnbay = ln(_bay = bay[j]); _nulny = nulny[j]

			// indexes in lambda[j,] after which switch from power or asymptotic-argument evaluation to asymptotic-parameter; min value is 0
			// Spuriously precise thresholds from benchmarking against mpmath's hyperu()
			_iy =  absb / (muj = _mu[j])
      _ix = editmissing(round(sqrt( (xhi[j]? _bax * 2.1846063 : (_bax * .01593272) ^ -1.1629507) * _iy - 1) / 1.6c782f6914edbX-003 /*h*/), 0)
			_iy = editmissing(round(sqrt( (yhi[j]? _bay * 2.1846063 : (_bay * .01593272) ^ -1.1629507) * _iy - 1) / 1.6c782f6914edbX-003 /*h*/), 0)
      anypowerx = xhi[j]==0 & _ix; anyasymargx = xhi[j] & _ix; anyasymparamx = _ix < $lnStickyFeller_N
      anypowery = yhi[j]==0 & _iy; anyasymargy = yhi[j] & _iy; anyasymparamy = _iy < $lnStickyFeller_N

			if (anypowerx) {  // prep for power series
				denompx1 = lnbax :- tp1
				denompx2 = lnbax :- tp2
			}
			if (anypowery) {
				denompy1 = lnbay :- tp1
				denompy2 = lnbay :- tp2
			}

			if (anyasymparamx)
				lnUAsymParamx = lnUAsymParam(alpha[|j,_ix+1\j,.|], beta, _bax, lnbax)
			if (anyasymparamy)
				lnUAsymParamy = lnUAsymParam(alpha[|j,_iy+1\j,.|], beta, _bay, lnbay)

			for (i=$lnStickyFeller_N;i;i--) {
        asymparamx = i > _ix; powerx = anypowerx & asymparamx==0; asymargx = anyasymargx & asymparamx==0
        asymparamy = i > _iy; powery = anypowery & asymparamy==0; asymargy = anyasymargy & asymparamy==0 

				_alpha = alpha[j,i]; lngammaalpha1beta = lngamma(_alpha1beta = _alpha + 1 - beta)
				lnC = _lnC + (lngammaalpha1beta - lngamma(_alpha))

				if (powerx | powery) {
					numerp1 = ln((_alpha      - 1) :+ lnStickyFeller_m)  // numerators in power series X::X+200.5 faster?
					numerp2 = ln((_alpha1beta - 1) :+ lnStickyFeller_m)
				}

				if (!(powerx & powery))
					lnM = lngammaalpha1beta - lngammanegnu

				if (asymargx | asymargy) {
					Ream1 = Re(_alpha) - 1; absa2 = abs(_alpha); absa2 = absa2*absa2; Ima2 = Im(_alpha); Ima2 = Ima2*Ima2
					ta = Ream1 - nu; ta2 = ta*ta + Ima2; ta3 = -Ream1 - Ream1 - 1 + absa2
					qa = 2*(Ream1+ta); qb = ta3+ta2+4*Ream1*ta; qc = 2*(Ream1*ta2+ta*ta3); qd = ta3*ta2  // coefficients of quartic equation governing when norm of ratio of successive terms is 1
					smaxx = asymargx? max((0, asdfQuarticRoot(qa, qb-_bax*_bax, qc, qd))) : -1 // power series index about where terms permanently start to increase
					smaxy = asymargy? max((0, asdfQuarticRoot(qa, qb-_bay*_bay, qc, qd))) : -1
					s = 0 :: max((smaxx, smaxy))
					numera = ln(s :+ (_alpha - 1)) + ln(s :+ (_alpha - beta)) - ln(s)  // numerator in large-argument asymptotic series; actually includes part for s! denominator too
				}

				retval = -ln(muj * lnStickyFeller_u2[i] / mu - nu * exp(lnC))
				if (_nulnx < .)
					if (powerx)
						retval = retval + hypergeomdiff(numerp1 + denompx1, lnC - _nulnx, numerp2 + denompx2)
					else if (asymparamx)
						retval = retval + lnUAsymParamx[i - _ix] + lnM
					else {  // asymptotic large-argument series
						(terms = (smaxx < smaxy? numera[|.\smaxx+1|] : numera) :- (1.921fb54442d18X+001i /*pi i*/ + lnbax)) [1] =  - _alpha * lnbax
						retval = retval + asdfLogSumExp(runningsum(terms)) + lnM
					}
				if (_nulny < .)
					if (powery)
						retval = retval + hypergeomdiff(numerp1 + denompy1, lnC - _nulny, numerp2 + denompy2)
					else if (asymparamy)
						retval = retval + lnUAsymParamy[i - _iy] + lnM
					else {  // asymptotic large-argument series
						(terms = (smaxy < smaxx? numera[|.\smaxy+1|] : numera) :- (1.921fb54442d18X+001i /*pi i*/ + lnbay)) [1] =  - _alpha * lnbay
						retval = retval + asdfLogSumExp(runningsum(terms)) + lnM
					}
				alpha[j,i] = retval					 
			}
		}

	} else {  // squared Bessel/b=0 case. Would need a different approach if b could vary by observation

		complex rowvector termx, lntermx, termy, lntermy, _lambda, _term, _lnterm, zerolimit; real scalar pidivsinnupi, gammanegnu2, c; real colvector ix, iy

		pidivsinnupi = pi() / sin(nu*pi())
		gammanegnu2 = gamma(-nu); gammanegnu2 = gammanegnu2 * gammanegnu2
		zerolimit = J(1, $lnStickyFeller_N, lngamma(-nu) - ln(2) * beta)
		c =                                            2 * ln(2) * beta 

		iy = (20.25 * a) :/ _mu  // 20.25 = 9^2/4 where 9 is the abs(z) cut-off between power series & asymptotic in BesselK(nu,z)
		ix = editmissing(round( sqrt(iy :/ x :- 1) / 1.6c782f6914edbX-003 /*h*/), 0)  // index of parabola points where BesselK argument passes 9
		iy = editmissing(round( sqrt(iy :/ y :- 1) / 1.6c782f6914edbX-003 /*h*/), 0)

		alpha = J(rows(x), $lnStickyFeller_N, C(.))
		for (j=rows(alpha);j;j--) {
			_term = sqrt(_lambda = _mu[j] * lnStickyFeller_u2) * (2 / sqrt(a)); _lnterm = ln(_term)  // precompute more of this
			termx = _term * sqrt(x[j]); lntermx = _lnterm :+ .5*lnx[j]
			termy = _term * sqrt(y[j]); lntermy = _lnterm :+ .5*lny[j]

			_ix = ix[j]; _iy = iy[j]
			alpha[j,] = c :+ ((x[j]? (_ix <  1                ? lnBesselKAsym (nu, termx, lntermx) : 
			                         (_ix >= $lnStickyFeller_N? lnBesselKPower(nu, lntermx)        :
			                                                    lnBesselKPower(nu, lntermx[|.\_ix|]), lnBesselKAsym(nu, termx[|_ix+1\.|], lntermx[|_ix+1\.|]))) - nu * lntermx :
																zerolimit) +
												(y[j]? (_iy < 1                 ? lnBesselKAsym (nu, termy, lntermy) : 
			                         (_iy >= $lnStickyFeller_N? lnBesselKPower(nu, lntermy)        : 
			                                                    lnBesselKPower(nu, lntermy[|.\_iy|]), lnBesselKAsym(nu, termy[|_iy+1\.|], lntermy[|_iy+1\.|]))) - nu * lntermy :
																zerolimit) - 
													 ln(_lambda * (gammanegnu2 / mu) - (_lambda / a) :^ -nu * pidivsinnupi))
		}
	}

	return (asdfLogSumExpRow((lnp0 , lnm + ln(_mu :* rowsum(asdfReexp(lnStickyFeller_tlambdalnu :+ alpha))))))
}


/// sum two hypergeometric series whose first term is 1 and is implicit. Multiply second by a coefficient, then subtract it from first. Avoid overflow. All inputs and output in logs
// *** overwrites x1, x2 for efficiency
complex scalar clsStickyFeller::hypergeomdiff(complex colvector x1, complex scalar c2, complex colvector x2) {
	real scalar shift
	                    _quadrunningsum(x1, x1)
	x2[1] = x2[1] + c2; _quadrunningsum(x2, x2)
	shift = ln(maxdouble()/rows(x1)) :- max((0, Re(c2), max(Re(x1)), max(Re(x2))))
	return (ln(exp(shift) - exp(c2 + shift) + quadcolsum(exp(x1 :+ shift) - exp(x2 :+ shift))) - shift)
}


// log Tricomi U for large a, dlmf.nist.gov/13.8.E11
// fails for entries of a whose real part is a negative integer because then lngamma(a) wrongly returns missing (as of 9/29/20)
complex rowvector clsStickyFeller::lnUAsymParam(complex rowvector a, real scalar b, real scalar z, real scalar lnz) {
  real colvector c, bzlnUAsymParam_B; real scalar k, m, t1, t2, t3, zs; complex rowvector term, lnterm, lnsqrtzdiva, lna; complex matrix mlna

	zs = k = 1
	c = (lnUAsymParam_q[1] = z/12 - b*.5) \ 1  // c_1 \ c_0 -- store c values backward
	bzlnUAsymParam_B = (-b) * lnUAsymParam_B - z * lnUAsymParam_kB

	for (m=1; m<=$lnUAsymParam_M-1; ) {
		t1 = bzlnUAsymParam_B[|.\++k|]
		c =  c          ' t1                  / k \ c  // add two more terms to c
		t1 = c[|.\k++|] ' t1
		c = (t1 + c[k] * bzlnUAsymParam_B[k]) / k \ c

		t1 = *lnUAsymParam_PT[m++] :* (zs = zs[1]*z \ zs)  // since c is stored backward, run these sums m to 0
		t3 = c[|.\m+1|] :/ (t2 = exp(lngamma(m+1-b :: .5-b)))  // exp(lngamma()) faster than gamma()!; even last entry is >0 if nu<0 so lngamma() will return right answer without converting argument to complex
		lnUAsymParam_p[m] = t2[2] * cross(t1, t3[|2\m+1|])
		lnUAsymParam_q[m] = t2[1] * cross(t1, t3[|.\m  |])
	}

	lna = ln(a)  // next 4 lines can be substantially precomputed
	lnsqrtzdiva = 0.5 * (lnz :- lna)
	term = sqrt(a * z); term = C(term + term); lnterm = ln(term)
	mlna = (0::$lnUAsymParam_M-1) * lna  // ln a^m
	return ((1.62e42fefa39efX-001 /*ln(2)*/ + z*.5) :+ ((1 - b) * lnsqrtzdiva - lngamma(a) + asdfLogSumExp(lnBesselKAsym(b-1, term, lnterm) :+ asdfLogSumExp(ln(C(lnUAsymParam_p)) :- mlna) \ 
	                                                                                                       lnBesselKAsym(b  , term, lnterm) :+ asdfLogSumExp(ln(C(lnUAsymParam_q)) :- mlna) :+ lnsqrtzdiva)))
}


// Compute log Bessel K with power series formula and complex argument z
// series length from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
complex rowvector clsStickyFeller::lnBesselKPower(real scalar nu, complex rowvector lnz) {
  real colvector InuDenom, InegnuDenom; complex matrix twomlnhalfz; complex rowvector nulnhalfz, lnhalfz

  (InuDenom    = ln(( nu::$lnBesselKPower_M.5+nu) :* lnBesselKPower_0toM)) [1] = lngamma(1+nu)
  (InegnuDenom = ln((-nu::$lnBesselKPower_M.5-nu) :* lnBesselKPower_0toM)) [1] = lngamma(1-nu)

	lnhalfz = lnz :- 1.62e42fefa39efX-001 /*ln .5*/
	nulnhalfz = nu * lnhalfz
	twomlnhalfz = lnBesselKPower_2Mto0 * lnhalfz
	return ((1.ce6bb25aa1315X-002 /*ln pi/2*/ - ln(C(sin(nu*pi())))) :+ asdfLogSumExp(asdfLogSumExp((twomlnhalfz :- nulnhalfz) :- runningsum(InegnuDenom)) \ 
	                                                                                  asdfLogSumExp((twomlnhalfz :+ nulnhalfz) :- runningsum(   InuDenom)) :+ 1.921fb54442d18X+001i /*pi i*/))
}


// Compute log Bessel K with asymptotic formula for small order (here, -1<nu<1 ASSUMED) and complex argument z
// max series length index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// This code only takes first 8 terms of series, which they recommend for abs(z) > 50.
complex rowvector clsStickyFeller::lnBesselKAsym(real scalar nu, complex rowvector z, complex rowvector lnz) {
  real scalar twonu, i; complex matrix terms; complex rowvector retval

  if (nu == -.5 | nu == .5)   // only first term is non-zero
		return (1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :- (.5*lnz + z))

	twonu = nu + nu
	(terms = J(1, i=cols(z), ln(C(twonu*twonu :- lnBesselKAsym_twokm1sq)) - lnBesselKAsym_k) :- lnz)[1,] = (-.5) * lnz - z  // log of first term, including outer multiplier exp(-z)
	retval = J(1,i,C(.))
	for (;i;i--)
		retval[i] = 1.ce6bb25aa1312X-003 + asdfLogSumExp(runningsum(terms[,i]))
	return (retval)	
}


complex rowvector asdfLogSumExp(complex matrix x) {
	real rowvector shift
	if (rows(x)==0) return(J(1,cols(x),0))
	if (rows(x)==1) return(x)
	shift = ln(maxdouble()/rows(x)) :- colmax(Re(x))
//	shift = shift - (shift:>0):*shift  // only downshift, to prevent overflow; shifting can prevent underflow & overflow but can also reduce precision if the shifter is much larger than entries
	return (ln(colsum(exp(x :+ shift))) - shift)
}


// Re(exp(z)) but a bit faster
real matrix asdfReexp(complex matrix z)	return (exp(Re(z)) :* cos(Im(z)))


// Find highest real root of a quartic equation using Euler's method, mathforum.org/dr.math/faq/faq.cubic.equations.html
real scalar asdfQuarticRoot(real scalar a, real scalar b, real scalar c, real scalar d) {
  real scalar e, f, fourg, h, i, j, a2, h2, h3; complex scalar alpha, z1, z2, p, q, r; complex colvector x

	    e = b - .375*(a2=a*a)  // substitute with y = x - b/4
	    f = c + (.125*a2 - .5*b)*a
	fourg = 4*d - a2*(.046875*a2 - .25*b) - a*c

	h = .5*e  // auxilliary cubic equation
	i = .0625*(e*e-fourg)
	j = -.015625*f*f

	p = (i-(h2=h*h)/3) / 3  // substite with z = y - h/3
	q = ((2*h2/27-i/3)*h+j) * .5

	alpha = (sqrt(C(q*q + p*p*p)) - q) ^ (1/3)
	z1 = alpha - p / alpha  // roots of p,q equation
	alpha = alpha * (-.5+1.bb67ae8584cabX-001i /*exp(2i*pi()/3)*/)
	z2 = alpha - p / alpha

	h3 = h/3
	p = sqrt(z1 - h3)  // square roots of roots of auxilliary cubic
	q = sqrt(z2 - h3)
	r = f/(p*q) * -.125

	x = p+q+r\p-q-r\-p+q-r\-p-q+r
	return (max(Re(select(x, Im(edittozerotol(x,1e-3)):==0)))-0.25*a)
}

mata mlib create lasdfStickyFeller, dir("`c(sysdir_plus)'l") replace
mata mlib add lasdfStickyFeller *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
----
* Goran test
mata S = clsStickyFeller(); p = exp(S.lnStickyFeller(t=rangen(.1,2,20)#J(20,1,1), x=J(400,1,1), y=J(20,1,1)#rangen(.1,2,20), a=1, b=1, nu=-.25, mu=5))
drop _all
getmata t y p, force double replace
twoway contour p y t, levels(500) clegend(off)
mata p :/ (exp(b/a*y) :* y:^nu :/ a)
