global M 200  // number of terms in power series expansion of M()/U() functions

mata
mata clear
mata set matastrict off
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
	complex rowvector lnUAsymParam(), lnUAsymParam2(), BesselKAsym(), lnBesselKAsym()
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


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007), Weideman (2010), parabolic contour method to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real colvector t, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu) {
  real colvector _mu, lnp0, negbt, lnattilde, lnm, lnx, lny, zerofill, zeroind, missingfill, lnJacobian; complex matrix lambda

	if (b > 0 | nu < -1 | nu > 0 | a<=0)
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

	_editmissing(lnm = (b / a) * y + nu * lny :- ln(a), -ln(mu))  // log speed measure; takes special value of 1/mu when y = 0

/*	log_epsilon = ln(1e-15)  // Logic extracted from Garrappa's mlf routine for Matlab
	log_eps = ln(epsilon(1))
	w = sqrt(log_eps / (log_eps - log_epsilon))  // half-width of integration range, needed to assure given precision
	N = ceil(-w * log_epsilon / (2*pi())) + 1    // half the number of integration points
	h = w :/ (N-.5)                              // width of bars in Riemann integral; do an even number of integration points since only integrating over half of parabola
  u = C(1, (.5..27.5) * 1.6c782f6914edbX-003)  // really 1 + ui */

  _mu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ :/ t  // different "mu" -- as in Garrappa (2015), p. 1364.
  lambda = _mu * lnStickyFeller_u2

  return (asdfLogSumExpRow((lnp0 , lnm + ln(_mu :* rowsum(asdfReexp(lnStickyFeller_tlambdalnu :+ lnStickyGreenTerm(lambda, x, y, a, b, nu, mu, lnx, lny)))))))
}

// Re(exp(z)) but a bit faster
real matrix asdfReexp(complex matrix z)	return (exp(Re(z)) :* cos(Im(z)))

// Compute sticky term of sticky Feller Green function for b<0. rows(lambda) must be max(rows(x),rows(y)) and rows(x) and rows(y) must each be that or 1
// accepts optional pre-computed ln(x), ln(y)
// MAY OVERWRITE lambda for efficiency
complex matrix clsStickyFeller::lnStickyGreenTerm(complex matrix lambda, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu, real colvector lnx, real colvector lny) {
	real colvector t1, t2, denomx1, denomx2, denomy1, denomy2, numer1, numer2, powerxj, poweryj
	complex matrix alpha
	real scalar _lnC, lnbax, lnbay, i, j, N, rowsx, rowsy, lngammanegnu, ba, bax, bay, powers, lnxnu, lnynu, onepnu
	real matrix powerx, powery
	complex scalar lnC, lnM, _alpha, _alphamnu, lngammaalphamnu
	complex rowvector lnUAsymParamx, lnUAsymParamy

	rowsx = rows(x); rowsy = rows(y)

	if (b < -1e-10) {
		powerxj = rowsum(powerx = abs(lambda :* (x/a)) :< 1.25e4f7b2737faX+005 /*53*ln(2)*/ )  // which entries to compute using with power series vs. asymptotic. Notably, twice this value will be fed into BesselK, where Zhang and Jin use cut-off of just 9 for asymptotic series.
		poweryj = rowsum(powery = abs(lambda :* (y/a)) :< 1.25e4f7b2737faX+005              )  // threshold from Johannson (2016), arXiv:1606.06977v2, p. 12; also from testing

		alpha = lambda / -b
		lngammanegnu = lngamma(-nu)
		_lnC = lngamma(onepnu = 1 + nu) - ln(-nu) - lngammanegnu - nu * ln(ba = -b/a)

		N = cols(lambda)

		if (any(powerxj) | any(poweryj)) {
			t1 = ln(onepnu :: $M.5 + nu) + lnStickyGreenTerm_lnm  // faster than ln(nu :+ m)
			t2 = ln(1 - nu :: $M.5 - nu) + lnStickyGreenTerm_lnm
		}

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

			if (any(powerxj[j]))
				lnUAsymParamx = lnUAsymParam(alpha[j,], onepnu, bax)
			if (any(poweryj[j]))
				lnUAsymParamy = lnUAsymParam(alpha[j,], onepnu, bay)

			for (i=N;i;i--) {
				_alpha = alpha[j,i]; _alphamnu = _alpha - nu; lngammaalphamnu = lngamma(_alphamnu)
				lnC = _lnC + (lngammaalphamnu - lngamma(_alpha))
				if (powers = powerx[j,i] + powery[j,i]) {
					numer1 = ln((_alpha    - 1) :+ lnStickyGreenTerm_m)  // X::X+200.5 faster?
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
  
	// squared Bessel/b=0 case. Would need a different approach if b could vary by observation

	complex rowvector termx, lntermx, termy, lntermy, _lambda, _term, _lnterm; real scalar gamma1pnudivgammanegnu, c, _x, _y, _lnx, _lny

	if (rowsx == 1) {
		_x = x
		_lnx = lnx
	}
	if (rowsy == 1) {
		_y = y
		_lny = lny
	}

	gamma1pnudivgammanegnu = gamma(1 + nu) / gamma(-nu)
	c = 2 * (1.62e42fefa39efX-001 /*ln(2)*/ - lngamma(-nu)) 
	
	for (j=rows(lambda);j;j--) {
		if (rowsx > 1) {
			_x = x[j]
			_lnx = lnx[j]
		}
		if (rowsy > 1) {
			_y = y[j]
			_lny = lny[j]
		}
		
		_term = sqrt(_lambda = lambda[j,]) * (2 / sqrt(a)); _lnterm = ln(_term)
		termx = _term * sqrt(_x); lntermx = _lnterm :+ .5*_lnx
		termy = _term * sqrt(_y); lntermy = _lnterm :+ .5*_lny
		lambda[j,] = (ln(BesselKAsym(nu, termx, lntermx)) + 
		              ln(BesselKAsym(nu, termy, lntermy)) -
		              nu * (lntermx + lntermy :- 1.62e42fefa39efX+000 /*ln(4)*/) - 
		              ln(_lambda / mu + (_lambda / a) :^ -nu * gamma1pnudivgammanegnu)) :+  c
	}
	return(lambda)
}

//  mata S = clsStickyFeller(); S.lnStickyGreenTerm(-5.72953504 - .236347915i,.001,.001,2,-2e-6,-.9,1,-6.90775528,-6.90775528)


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
  real colvector c, bzlnUAsymParam_B; real rowvector p, q; real scalar k, s, _q, t1, t2, t3, zs; complex rowvector term, lnterm, C1, C2, zdiva; complex matrix as

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
		t3 = c[|.\s+2|] :/ (t2 = exp(lngamma(s+2-b :: .5-b)))  // exp(lngamma()) faster than gamma()!; note that even last entry is >0 if nu<0 so lngamma() will return right answer without converting argument to complex
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


// log Tricomi U for large a, dlmf.nist.gov/13.8.E11
// fails for entries of a whose real part is a negative integer because then lngamma(a) wrongly returns missing (as of 9/29/20)
complex rowvector clsStickyFeller::lnUAsymParam2(complex rowvector a, real scalar b, real scalar z) {
  real colvector c, bzlnUAsymParam_B; real rowvector p, q; real scalar k, s, t1, t2, t3, zs; complex rowvector term, lnterm, lnsqrtzdiva, lna; complex matrix slna

	p = q = J(12,1,1)  // * should precompute
	zs = k = 1
	c = (z/12 - b*.5) \ 1  // c_1 \ c_0 -- store c values backward
	bzlnUAsymParam_B = (-b) * lnUAsymParam_B - z * lnUAsymParam_kB

	for (s=1; s<=11; s++) {
		t1 = bzlnUAsymParam_B[|.\++k|]
		c =  c          ' t1                  / k \ c  // add two more terms to c
		t1 = c[|.\k++|] ' t1
		c = (t1 + c[k] * bzlnUAsymParam_B[k]) / k \ c

		t1 = *lnUAsymParam_PT[s] :* (zs = zs[1]*z \ zs)  // since c is stored backward, run these sums s to 0
		t3 = c[|.\s+2|] :/ (t2 = exp(lngamma(s+2-b :: .5-b)))  // exp(lngamma()) faster than gamma()!; note that even last entry is >0 if nu<0 so lngamma() will return right answer without converting argument to complex
		p[s+1] = t2[2] * cross(t1, t3[|2\s+2|])
		q[s+1] = t2[1] * cross(t1, t3[|.\s+1|])
	}
	
	lna = ln(a)  // next 4 lines can be substantially precomputed
	lnsqrtzdiva = 0.5 * (ln(z) :- lna)
	term = sqrt(a * z); term = C(term + term); lnterm = ln(term)
	slna = J(1,cols(a),(0::11) * lna)  // ln a^s  // in general cols(a) is fixed (28)
	return ((1.62e42fefa39efX-001 /*ln(2)*/ + z*.5) :+ ((1 - b) * lnsqrtzdiva - lngamma(a) + asdfLogSumExp(lnBesselKAsym(b-1, term, lnterm) :+ asdfLogSumExp(ln(p) :- slna) \ 
	                                                                                                        lnBesselKAsym(b  , term, lnterm) :+ asdfLogSumExp(ln(q) :- slna) :+ lnsqrtzdiva)))
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

// Compute log Bessel K with asymptotic formula for small order (here, -1<nu<1 ASSUMED) and complex argument z
// max series length index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// *** This code only takes first 8 terms of series, which they recommend for abs(z) > 50.
complex rowvector clsStickyFeller::lnBesselKAsym(real scalar nu, complex rowvector z, complex rowvector lnz) {
  real scalar twonu, i; complex matrix terms; complex rowvector retval

  if (nu == -.5 | nu == .5)   // only first term is non-zero
		return (1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :- (.5*lnz + z))

	twonu = nu + nu
	(terms = J(1, cols(z), ln(twonu*twonu :- BesselKAsym_twokm1sq) :- BesselKAsym_k) :- lnz)[1,] = (-.5) * lnz - z  // log of first term, including outer multiplier exp(-z)
	if (cols(z)==1)
		return (    1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :+ asdfLogSumExp(runningsum(terms    )))

	retval = J(1, cols(z), C(.))  // note that this will always have same width--currently 27
	for (i=cols(z);i;i--)
		retval[i] = 1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :+ asdfLogSumExp(runningsum(terms[,i]))
	return (retval)	
}


complex rowvector asdfLogSumExp(complex matrix x) {
	real rowvector shift
	if (rows(x)==0) return(J(1,cols(x),0))
	shift = ln(maxdouble()/rows(x)) :- colmax(Re(x))
//	shift = shift - (shift:>0):*shift  // only downshift, to present overflow; shifting can prevent underflow & overflow but can also reduce precision if the shifter is much larger than entries
	return (ln(colsum(exp(x :+ shift))) - shift)
}


// Find highest root of a quartic equation using Euler's method, http://mathforum.org/dr.math/faq/faq.cubic.equations.html
complex colvector asdfQuarticRoots(real scalar a, real scalar b, real scalar c, real scalar d) {
  real scalar e, f, g, h, i, j, a2, h2; complex scalar alpha, z1, z2, p, q, r

	e = b - .375*(a2=a*a)  // substitute with y = x - b/4
	f = c + .125*a2*a - .5*a*b
	g = d - .01171875*a2*a2 + .0625*a2*b - .25*a*c

	h = 0.5*e  // auxilliary cubic equation
	i = .0625*(e*e-4*g)
	j = -.015625*f*f

	p = (i-(h2=h*h)/3) / 3  // substite with z = y - h/3
	q = (2*h2*h/27-i*h/3+j) / 2

	alpha = (sqrt(C(q*q + p*p*p)) - q) ^ (1/3)
	z1 = alpha - p / alpha  // roots of p,q equation
	alpha = alpha * exp(2i*pi()/3)
	z2 = alpha - p / alpha

	p = sqrt(z1 - h/3)  // square roots of roots of auxilliary cubic
	q = sqrt(z2 - h/3)
	r = f/(p*q) * -.125

	return ((p+q+r\p-q-r\-p+q-r\-p-q+r):-0.25*a)
}


// log Tricomi U for large z, dlmf.nist.gov/13.7.E3
// takes ln z instead of z
complex scalar clsStickyFeller::lnUAsymArg(complex scalar a, real scalar b, real scalar lnz) {
	real colvector s; complex colvector terms
	s = -1 :: max(Re(asdfQuarticRoots(2*(2*Re(a)-b+1), 
	                                  abs(a)^2+(Re(a)-b+1)^2+Im(a)^2+4*Re(a)*(Re(a)-b+1)-exp(2*lnz), 
								                    2*(Re(a)*((Re(a)-b+1)^2+Im(a)^2)+(Re(a)-b+1)*abs(a)^2),
								                    abs(a)^2*((Re(a)-b+1)^2+Im(a)^2))))

	(terms = ln(s :+ a) + ln(s :+ (a - b + 1)) - ln(s:+1) :- (1.921fb54442d18X+001i /*pi i*/ + lnz)) [1] =  - a * lnz
	return (asdfLogSumExp(runningsum(terms)))
}

mata mlib create lasdfStickyFeller, dir("`c(sysdir_plus)'l") replace
mata mlib add lasdfStickyFeller *(), dir("`c(sysdir_plus)'l")
mata mlib index
end

cap python set userpath "C:\Users\drood\AppData\Local\Programs\Python\Python38-32\lib\site-packages", perm
python
import mpmath as mp
from sfi import *
end

cap mata mata drop M()
cap mata mata drop lnUPower()
cap mata mata drop lnU2F0()

mata
function lnU2F0(alpha,beta,z) {
	m = 0::200
	retval = lngamma(alpha:+m)+lngamma((alpha-beta+1):+m)-lngamma(m:+1) -m*ln(C(-z))
 i = 1 + min((1::200) :/ (Re(retval[|.\200|]) :< Re(retval[|2\.|])))  // crude way to terminate series once it starts to increase
/*"(-alpha*ln(z)-lngamma(alpha)-lngamma(alpha-beta+1)):+retval"
 (-alpha*ln(z)-lngamma(alpha)-lngamma(alpha-beta+1)):+retval
"i"
i*/
	return (-alpha*ln(z)-lngamma(alpha)-lngamma(alpha-beta+1) + asdfLogSumExp(retval[|.\i|]))
}


function M(alpha,beta,z) {
	real colvector m
	m = 0::2000
	return (quadcolsum(exp(lngamma(alpha:+J(1,cols(alpha),m)) :- lngamma(beta:+m) :- lngamma(m:+1) :+ m :* ln(z) :- lngamma(alpha) :+ lngamma(beta))))
}
function lnUPower (alpha,beta,z)
	return(asdfquadLogSumExp(colsum(                lngamma(  1-beta) :-lngamma(alpha:+1:-beta):+ln(M(alpha,            beta,z))) \ 
	                         colsum((1-beta)*ln(z) +lngamma(C(beta-1)):-lngamma(alpha         ):+ln(M(alpha:+1:-beta,2:-beta,z)))))

mata
S = clsStickyFeller()
lambda=-1+1i;a=1.1;b=-1;nu=-.75;x=1.2
alpha = -lambda/b; beta=1+nu; z=-b*x/a

table = J(0,7,0)
for (i=1;i<=10;i=i+1/2)
	for (j=1;j<=10;j=j+1/2) {
		_alpha = alpha*exp(i-2)
		_z = z*exp(j-5)
i,j,_alpha,_z
		t = C(.)
		stata(`"cap python: Mata.store("t", complex(mp.ln(mp.hyperu(Mata.getAt("_alpha",0,0), 1+Mata.getAt("nu",0,0), Mata.getAt("_z",0,0)))))"')
		table = table \ _alpha, _z, lnUPower(_alpha,beta,_z), S.lnUAsymParam2(_alpha,beta,_z), S.lnUAsymArg(_alpha,beta,ln(_z)) , lnU2F0(_alpha,beta,_z), t
	}
table
Retable = Re(table)
best=abs(Retable[,3..6] :- Retable[,7]); best=rowsum((best:==rowmin(best)) :* (1..4))
end
getmata (alpha z lnU lnUAsymParam lnUAsymArg lnU2F0 hyperu) = Retable, force double replace
cap label drop functions
label define functions 1 "lnU" 2 "Param" 3 "Arg" 4 "2F0"
getmata best, force replace
label values best functions
replace alpha = -alpha
replace best = . if best==10
scatter alpha z if best==1, xscale(log) yscale(log) scheme(s1rcolor) mlab(best) || ///
	scatter alpha z if best==1,  mlab(best) || scatter alpha z if best==2,  mlab(best) || ///
	scatter alpha z if best==3,  mlab(best) || scatter alpha z if best==4,  mlab(best) legend(off)