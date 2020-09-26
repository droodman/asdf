global M 200  // number of terms in power series expansion of M()/U() functions

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

class clsStickyFeller {
	real colvector lnStickyGreenTerm_m, lnStickyGreenTerm_lnm, lnUAsymParam_B, lnUAsymParam_kB
	real colvector BesselKAsym_k, BesselKAsym_twokm1sq
	complex rowvector lnStickyFeller_tlambdalnu, lnStickyFeller_u2
	pointer(real colvector) colvector lnUAsymParam_PT

	real colvector lnStickyFeller()
	complex matrix lnStickyGreenTerm(), lnStickyGreen()
	complex scalar lnUAsymArg(), hypergeomdiff()
	complex rowvector lnUAsymParam(), BesselKAsym()
	void new()
}

void clsStickyFeller::new() {
	lnStickyGreenTerm_m = 1::$M; lnStickyGreenTerm_lnm = ln(lnStickyGreenTerm_m)

  complex rowvector u
	u = C(1, (.5..27.5) * 1.6c782f6914edbX-003 /*h*/)  // really (1 + ui)
	lnStickyFeller_u2 = u :* u
	lnStickyFeller_tlambdalnu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ * lnStickyFeller_u2 + ln(u) :- 1.16c0cdab51479X+001 /* + ln 2*2*h/2pi  tack this on here for effciency */

	BesselKAsym_k = 0::7  // 7 comes from using 1st 8 terms of series
	BesselKAsym_twokm1sq = 2*BesselKAsym_k:-1; BesselKAsym_twokm1sq = BesselKAsym_twokm1sq :* BesselKAsym_twokm1sq
	BesselKAsym_k = ln(8 * BesselKAsym_k)  // 8 comes from asymptotic formula

	// B_i / i!, i=1 to 20, B-i = Bernoulli numbers
	lnUAsymParam_B = -1.0000000000000X-001 \ 1.5555555555555X-004 \ 0 \ -1.6c16c16c16c17X-00a \ 0 \ 1.1566abc011566X-00f \ 0 \ -1.bbd779334ef0bX-015 \ 0 \ 1.66a8f2bf70ebeX-01a \ 0 \ -1.22805d644267fX-01f \ 0 \ 1.d6db2c4e09163X-025 \ 0 \ -1.7da4e1f79955cX-02a \ 0 \ 1.355871d652e9eX-02f \ 0 \ -1.f57d968caacf1X-035 \ 0 \ 1.967e1f09c376fX-03a \ 0 \ -1.39544646858125233407076862640635497639176366691109952378E-19 \ 0 \ 3.5347070396294674716932299778037992147245945647149203595069E-21 \ 0 \ -8.95351742703754685040261131811274105162713924278496251644E-23 \ 0 \ 2.2679524523376830603109507388681660632203543297441712610647E-24 // (-1/2\ 1/6\ 0 \-1/30\ 0 \1/42\ 0 \-1/30\ 0 \5/66\ 0 \-691/2730\ 0 \7/6\ 0 \-3617/510\ 0 \43867/798\ 0 \-174611/330 \ 0 \ 854513/138) :/ factorial(1::22)
	lnUAsymParam_kB = lnUAsymParam_B[|2\.|] :* (1::rows(lnUAsymParam_B)-1) \ 0

	// Pascal's triangle
	lnUAsymParam_PT = &(1\1)
	real matrix PT; real scalar s
	PT = comb(2..20,0::20)
	for (s=1;s<=20-2;s++)
		lnUAsymParam_PT = lnUAsymParam_PT \ &(PT[|.,s\s+2,s|])
//	lnUAsymParam_PT = &(1\1) \ &(1\2\1) \ &(1\3\3\1) \ &(1\4\6\4\1) \ &(1\5\10\10\5\1) \ &(1\6\15\20\15\6\1) \ &(1\7\21\35\35\21\7\1) \ &(1\8\28\56\70\56\28\8\1) \ &(1\9\36\84\126\126\84\36\9\1) \ &(1\10\45\120\210\252\210\120\45\10\1) \ &(1\11\55\165\330\462\462\330\165\55\11\1)
}


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007) parabolic contour method just to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real colvector t, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu) {
  real colvector _mu, lnp0, negbt, lnattilde, lnm, lnx, lny, zerofill, zeroind, missingfill, lnJacobian; complex matrix lambda

	if (b >= 0 | nu < -1 | nu > 0 | a<=0)
		return (.)

	zerofill = J(length(zeroind = selectindex(y:==0)), 1, 0)
	missingfill = J(length(zeroind), 1, .)
	lnx = ln(x); lny = ln(y)

	lnattilde = ln(a / (-b) * expm1(negbt = (-b) * t))
	lnJacobian =  negbt - lnattilde
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

	_editmissing(lnm = (b / a) * y + nu * lny :- ln(a), -ln(mu))  // log speed measure; takes special value of 1/mu when y = 0

/*	log_epsilon = ln(1e-15)  // Logic extracted from Garrappa's mlf routine for Matlab
	log_eps = ln(epsilon(1))
	w = sqrt(log_eps / (log_eps - log_epsilon))  // half-width of integration range, needed to assure given precision
	N = ceil(-w * log_epsilon / (2*pi())) + 1    // half the number of integration points
	h = w :/ (N-.5)                              // width of bars in Riemann integral; do an even number of integration points since only integrating over half of parabola
  u = C(1, (.5..27.5) * 1.6c782f6914edbX-003)  // really 1 + ui */

  _mu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ :/ t  // different "mu" -- as in Weideman and Trefethen eq 1.9
  lambda = _mu * lnStickyFeller_u2
  return (asdfLogSumExpRow((lnp0 , lnm + ln(_mu :* rowsum(asdfReexp(lnStickyFeller_tlambdalnu :+ lnStickyGreenTerm(lambda, x, y, a, b, nu, mu, lnx, lny)))))))
}

// Re(exp(z)) but a little faster
real matrix asdfReexp(complex matrix z)	return (exp(Re(z)) :* cos(Im(z)))

// Compute sticky term of sticky Feller Green function for b<0. rows(lambda) must be max(rows(x),rows(y)) and rows(x) and rows(y) must each be that or 1
// accepts optional pre-computed ln(x), ln(y)
complex matrix clsStickyFeller::lnStickyGreenTerm(complex matrix lambda, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu, | real colvector lnx, real colvector lny) {
	real colvector t1, t2, denomx1, denomx2, denomy1, denomy2, numer1, numer2, powerxj, poweryj
	complex matrix alpha
	real scalar _lnC, lnbax, lnbay, i, j, N, rowsx, rowsy, lngammanegnu, ba, bax, bay, powers, lnxnu, lnynu, onepnu
	real matrix powerx, powery
	complex scalar lnC, lnM, _alpha, _alphamnu, lngammaalphamnu
	complex rowvector lnUAsymParamx, lnUAsymParamy

	powerxj = rowsum(powerx = abs(lambda :* (x/a)) :< 30)  // which entries to compute using with power series vs. asymptotic. Notably, twice this value will be fed into BesselK, where Zhang and Jin use cut-off of just 9 for asymptotic series. Should this threshold go as low as 4.5?
	poweryj = rowsum(powery = abs(lambda :* (y/a)) :< 30)
	alpha = lambda / -b
	lngammanegnu = lngamma(-nu)
	_lnC = lngamma(onepnu = 1 + nu) - ln(-nu) - lngammanegnu - nu * ln(ba = -b/a)

	N = cols(lambda)
	rowsx = rows(x); rowsy = rows(y)
	if (rows(lnx) == 0) lnx = ln(x)
	if (rows(lny) == 0) lny = ln(y)

	t1 = ln(onepnu :: $M.5 + nu) + lnStickyGreenTerm_lnm  // faster than ln(nu :+ lnStickyGreenTerm_m)
	t2 = ln(1 - nu :: $M.5 - nu) + lnStickyGreenTerm_lnm

	if (rowsx == 1) {
		lnxnu = nu * lnx
		bax = ba * x
		if (any(powerxj)) {
			lnbax = ln(bax)
			denomx1 = lnbax :- t1
			denomx2 = lnbax :- t2
		}
	}
	
	if (rowsy == 1) {
		lnynu = nu * lny
		bay = ba * y
		if (any(poweryj)) {
			lnbay = ln(bay)
			denomy1 = lnbay :- t1
			denomy2 = lnbay :- t2
		}
	}
	
	for (j=rows(lambda);j;j--) {
		if (rowsx > 1) {
			lnxnu = nu * lnx[j]
			bax = ba * x[j]
			if (powerxj[j]) {
				lnbax = ln(bax)
				denomx1 = lnbax :- t1
				denomx2 = lnbax :- t2
			}
		}
		if (rowsy > 1) {
			lnynu = nu * lny[j]
			bay = ba * y[j]
			if (poweryj[j]) {
				lnbay = ln(bay)
				denomy1 = lnbay :- t1
				denomy2 = lnbay :- t2
			}
		}

		if (all(powerxj) == 0)  // possibly suboptimal: compute asymptotic U~ for all x or y if any require it
			lnUAsymParamx = lnUAsymParam(alpha[j,], onepnu, bax)
		if (all(poweryj) == 0)
			lnUAsymParamy = lnUAsymParam(alpha[j,], onepnu, bay)

		for (i=N;i;i--) {
			_alpha = alpha[j,i]; _alphamnu = _alpha - nu; lngammaalphamnu = lngamma(_alphamnu)
			lnC = _lnC + (lngammaalphamnu - lngamma(_alpha))
			if (powers = powerx[j,i] + powery[j,i]) {
				numer1 = ln((_alpha    - 1) :+ lnStickyGreenTerm_m)
				numer2 = ln((_alphamnu - 1) :+ lnStickyGreenTerm_m)
			}
			if (powers < 2)
				lnM = lngammaalphamnu - lngammanegnu
			alpha[j,i] = (lnxnu == .? 0 : (powerx[j,i]? hypergeomdiff(numer1 + denomx1, lnC - lnxnu, numer2 + denomx2) : lnUAsymParamx[i] + lnM)) + 
			             (lnynu == .? 0 : (powery[j,i]? hypergeomdiff(numer1 + denomy1, lnC - lnynu, numer2 + denomy2) : lnUAsymParamy[i] + lnM)) -
			               ln(lambda[j,i] / mu - nu * exp(lnC))
		}
	}
	return(alpha)
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
complex rowvector clsStickyFeller::lnUAsymParam(complex rowvector a, real scalar b, real scalar z) {
  real colvector c, bzlnUAsymParam_B; real rowvector p, q; real scalar k, s, _q, t1, t2, t3, zs; complex rowvector retval, term, lnterm, C1, C2, zdiva; complex matrix as

	p = q = J(1,11,0)
	zs = k = 1
	c = (_q = z/12 - b*.5) \ 1  // c_1 \ c_0 -- store c values backward
	bzlnUAsymParam_B = (-b) * lnUAsymParam_B - z * lnUAsymParam_kB

	for (s=1; s<=11; s++) {
		t1 = bzlnUAsymParam_B[|.\++k|]
		c =  c          ' t1                  / k \ c  // add two more terms to c
		t1 = c[|.\k++|] ' t1
		c = (t1 + c[k] * bzlnUAsymParam_B[k]) / k \ c

		t1 = *lnUAsymParam_PT[s] :* (zs = zs[1]*z \ zs)  // since c is stored backward, run these sums s to 0
		t3 = c[|.\s+2|] :/ (t2 = exp(lngamma(s+2-b :: .5-b)))  // exp(lngamma()) faster than gamma()!
		p[s] = t2[2] * cross(t1, t3[|2\s+2|])
		q[s] = t2[1] * cross(t1, t3[|.\s+1|])
	}
	
	zdiva = z :/ a
	term = sqrt(a * z); term = C(term + term)
	lnterm = ln(term)
	as = exp((-1::-11) * ln(a))  // a^-s
	C1 =                     BesselKAsym(b-1, term, lnterm)
	C2 = _q * sqrt(zdiva) :* BesselKAsym(b  , term, lnterm)
	return ((1.62e42fefa39efX-001 /*ln(2)*/  + z*.5) :+ (1 - b)*.5 * ln(zdiva) + ln(C1 + C1 :* p * as + C2 + C2 :* q * as) :- lngamma(a))
}


// log Tricomi U for large z, dlmf.nist.gov/13.7.E3
// takes ln z instead of z
complex scalar clsStickyFeller::lnUAsymArg(complex scalar a, real scalar b, real scalar lnz) {
	real colvector s; complex colvector terms
	s = -1::8  // 10 terms is made-up; can be precomputed; or need to iterate and check for explosion?
	(terms = ln(s :+ a) + ln(s :+ (a - b + 1)) - editmissing(ln(s),0) - (s:+1) * (1.921fb54442d18X+001i /*pi i*/ + lnz))[1] =  - a * lnz
	return (asdfLogSumExp(runningsum(terms)))
}


// Compute Bessel K with asymptotic formula for small order (here, -1<nu<1 ASSUMED) and complex argument z
// max series length index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// *** This code only takes first 8 terms of series, which they recommend for abs(z) > 50.
complex rowvector clsStickyFeller::BesselKAsym(real scalar nu, complex rowvector z, complex rowvector lnz) {
  real scalar twonu, i; complex matrix terms; complex rowvector retval

  if (nu == -.5 | nu == .5)   // only first term is non-zero
		return (1.40d931ff62705X+000 /*sqrt(pi/2)*/ :/ (sqrt(z) :* exp(z)))

	twonu = nu + nu
	(terms = J(1, cols(z), ln(twonu*twonu :- BesselKAsym_twokm1sq) :- BesselKAsym_k) :- lnz)[1,] = (-.5) * lnz - z  // log of first term, including outer multiplier exp(-z)
	if (cols(z)==1)
		return (    1.40d931ff62705X+000 /*sqrt(pi/2)*/ * colsum(exp(runningsum(terms   , 1))))

	retval = J(1, cols(z), C(.))
	for (i=cols(z);i;i--)
		retval[i] = 1.40d931ff62705X+000 /*sqrt(pi/2)*/ * colsum(exp(runningsum(terms[,i], 1)))
	return (retval)	
}

mata mlib create lasdfStickyFeller, dir("`c(sysdir_plus)'l") replace
mata mlib add lasdfStickyFeller *(), dir("`c(sysdir_plus)'l")
mata mlib index
end