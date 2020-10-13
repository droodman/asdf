global UPower_M 60  // max term index in power series expansion of U()
global BesselKPower_M 22  // max term index in power series expansion of BesselK()
global BesselKAsym_M 7  // max term index in asymptotic approximation of BesselK()
global UAsymParam_M 12 // max term index in large-parameter approximation of U()
global StickyFeller_N 28  // number of points sampled on parabola of integration

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

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

// Compute log Bessel K with asymptotic formula for small order (here, -1<nu<1 ASSUMED) and complex argument z
// max series length index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// This code only takes first 8 terms of series, which they recommend for abs(z) > 50.
complex rowvector clsBesselKAsym::lnK(real scalar nu, complex rowvector z, complex rowvector lnz) {
  real scalar twonu, i; complex matrix terms; complex rowvector retval

  if (nu == -.5 | nu == .5)   // only first term is non-zero
		return (1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :- (.5*lnz + z))

	twonu = nu + nu
	(terms = J(1, i=cols(z), ln(C(twonu*twonu :- twokm1sq)) - k) :- lnz)[1,] = (-.5) * lnz - z  // log of first term, including outer multiplier exp(-z)
	retval = J(1,i,C(.))
	for (;i;i--)
		retval[i] = 1.ce6bb25aa1312X-003 + asdfLogSumExp(quadrunningsum(terms[,i]))
	return (retval)	
}

class clsUPower {
  pointer(complex rowvector) scalar palpha
	real scalar beta, paramDirty, lnz, i
  real colvector m, lnm, mlnz, denom1, denom2, mbeta
  complex colvector terms1, terms2
  complex matrix coefs1, coefs2
  
  void new(), setalpha(), setbeta(), setz()
  complex rowvector lnU()
}

void clsUPower::new() {
	lnm  = ln(m = 0::$UPower_M)
  coefs1 = coefs2 = J(`=$UPower_M+1', $StickyFeller_N, C(.))
}

void clsUPower::setalpha(complex rowvector _alpha) {
	if (palpha==NULL)
  	palpha = &_alpha
  else if (cols(*palpha) != cols(_alpha) | (*palpha)[1] != _alpha[1])
  	palpha = &_alpha
  paramDirty = 1
}

void clsUPower::setbeta(real scalar _beta) {
  if (beta != _beta) {
    beta = _beta
    paramDirty = 1
    denom1 = quadrunningsum(ln(        beta - 1 :: $UPower_M.5 + beta - 1) + lnm)  // first entry will be 0 since lnm[1]=.
    denom2 = quadrunningsum(ln(mbeta = 1 - beta :: $UPower_M.5 + 1 - beta) + lnm)
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
complex rowvector clsUPower::lnU(real scalar maxi /*complex rowvector lnC*/) {
  complex scalar _alpha; real colvector shift; complex rowvector S1, S2

  if (maxi > $StickyFeller_N) maxi = $StickyFeller_N
  for (i=paramDirty? 1 : i+1; i <= maxi; i++) {
    _alpha = (*palpha)[i]
    (terms1 = ln((_alpha -    1) :+ m)) [1] = 0       // for making log Pochhammer symbols (alpha)_m and (alpha+1-beta)_m; edit m=0 entry to ln 1
    (terms2 = ln((_alpha - beta) :+ m)) [1] = lngamma(beta)-lngamma(2-beta)+lngamma(_alpha+1-beta)-lngamma(_alpha)  // edit m=0 entry to multiplier on second series
    coefs1[,i] = quadrunningsum(terms1) - denom1
    coefs2[,i] = quadrunningsum(terms2) - denom2
  }
  paramDirty = 0
  S1 = asdfLogSumExp(coefs1 :+        mlnz)
  S2 = asdfLogSumExp(coefs2 :+ mbeta * lnz)
  shift = 1.6232bdd7abcd2X+009 /*ln(maxdouble()/2)*/ :- colmax(Re(S1) \ Re(S2))
  return (ln(exp(S1 + shift) - exp(S2 + shift)) - shift)  // return power series sum in logs
}


class clsUAsymParam {
  pointer(complex rowvector) scalar palpha
	real scalar beta, betaDirty, zDirty, z, lnz
	real colvector B, kB, zmInd
	real rowvector mlnz
	pointer(real rowvector) colvector pPT
	complex rowvector lnalpha, lngammaalpha, lnSp, lnSq
	complex matrix lnp, lnq, C2, mlnalpha
  class clsBesselKAsym scalar SBessel

  void new(), setalpha(), setbeta(), setz()
  complex rowvector lnU()
}

void clsUAsymParam::new() {
	real matrix PT; real scalar s

	// B_i / i!, i=1 to 99, B_i = Bernoulli numbers
	B = -.5 \ 0.0833333333333333 \ 0 \ -.00138888888888889 \ 0 \ 3.30687830687831e-5 \ 0 \ -8.26719576719577e-7 \ 0 \ 2.08767569878681e-8 \ 0 \ -5.28419013868749e-10 \ 0 \ 1.33825365306847e-11 \ 0 \ -3.38968029632258e-13 \ 0 \ 8.58606205627785e-15 \ 0 \ -2.17486869855806e-16 \ 0 \ 5.50900282836023e-18 \ 0 \ -1.39544646858125e-19 \ 0 \ 3.53470703962947e-21 \ 0 \ -8.95351742703755e-23 \ 0 \ 2.26795245233768e-24 \ 0 \ -5.7447906688722e-26 \ 0 \ 1.45517247561486e-27 \ 0 \ -3.68599494066531e-29 \ 0 \ 9.33673425709505e-31 \ 0 \ -2.36502241570063e-32 \ 0 \ 5.99067176248213e-34 \ 0 \ -1.51745488446829e-35 \ 0 \ 3.84375812545419e-37 \ 0 \ -9.73635307264669e-39 \ 0 \ 2.46624704420068e-40 \ 0 \ -6.24707674182074e-42 \ 0 \ 1.58240302446449e-43 \ 0 \ -4.00827368594894e-45 \ 0 \ 1.01530758555696e-46 \ 0 \ -2.57180415824187e-48 \ 0 \ 6.51445603523381e-50 \ 0 \ -1.65013099068965e-51 \ 0 \ 4.17983062853948e-53 \ 0 \ -1.05876346677029e-54 \ 0 \ 2.68187919126077e-56 \ 0 \ -6.79327935110742e-58 \ 0 \ 1.72075776166814e-59 \ 0 \ -4.35873032934889e-61 \ 0 \ 1.10407929036847e-62 \ 0 \ -2.79666551337813e-64 \ 0 \ 7.08403650167947e-66 \ 0 \ -1.79440740828922e-67 \ 0 \ 4.5452870636111e-69 \ 0 \ -1.15133466319821e-70 \ 0 \ 2.91636477109236e-72 \ 0 \ -7.38723826349734e-74 \ 0 \ 1.8712093117638e-75 \ 0 \ -4.7398285577618e-77 \ 0 \ 1.20061259933545e-78 \ 0
	B = -B  // in lmf.nist.gov/13.8.E16, moving sum to RHS negates terms 
  kB = B[|2\.|] :* (1::rows(B)-1) \ 0

	// Pascal's triangle
	pPT = J($UAsymParam_M+1, 1, NULL)
	PT = comb(0::`=$UAsymParam_M+1', 0..`=$UAsymParam_M+1')
	for (s=1; s<=`=$UAsymParam_M+1'; s++)
		pPT[s] = &(PT[|s,.\s,s|])
	
	lnp = lnq = J(`=$UAsymParam_M+1', `=3*$UAsymParam_M+2', C(.))
	lnp[1,1] = 0; lnq[1,2] = 1.3e116bcd39e7dX+001 /* ln 1/12 */  // first-order term of q_0(z) is always z/12

	zmInd = (1::`=(2*$UAsymParam_M+2)^2') + (`=2*(2*$UAsymParam_M+2)^2+2*$UAsymParam_M+2') :- `=3*(2*$UAsymParam_M+2)'*(1::`=2*$UAsymParam_M+2')#J(`=2*$UAsymParam_M+2',1,1)  // used to premultiply c_m polys by z^m and flip order
	C2 = J(`=(2*$UAsymParam_M+2)*(2*(2*$UAsymParam_M+2)-1)', 1, 0)
	
	SBessel = clsBesselKAsym() 
}

void clsUAsymParam::setalpha(complex rowvector _alpha) {  // maybe can be dispensed with through "extends"
	if (palpha==NULL) {
  	palpha = &_alpha
		lnalpha = ln(_alpha)
		mlnalpha = (0 :: $UAsymParam_M) * lnalpha
		lngammaalpha = lngamma(_alpha)
	} else if (cols(*palpha) != cols(_alpha) | (*palpha)[1] != _alpha[1]) {
  	palpha = &_alpha
		lnalpha = ln(_alpha)
		mlnalpha = (0 :: $UAsymParam_M) * lnalpha
		lngammaalpha = lngamma(_alpha)
	}
}

void clsUAsymParam::setbeta(real scalar _beta) {
  if (beta != _beta) {
    beta = _beta
    betaDirty = 1
  }
}

void clsUAsymParam::setz(real scalar _z, real scalar _lnz) {
	if (lnz != _lnz) {
    z    = _z
    lnz  = _lnz
		zDirty = 1
  }
}

// log Tricomi U for large a, dlmf.nist.gov/13.8.E11
// fails for entries of a whose real part is a negative integer because then lngamma(a) wrongly returns missing (as of 9/29/20)
complex rowvector clsUAsymParam::lnU() {
  real colvector betaB; real scalar m, t1, t2; complex rowvector term, lnterm, lnsqrtzdivalpha; real matrix C, t3

	if (betaDirty) {
		C = 1 // c_0 -- store c polynomials bottom-to-top
		betaB = beta * B  // constant terms of binomials inside dlmf.nist.gov/13.8.E16; 1st-order terms already in kB

		for (m=1; m<=`=2*$UAsymParam_M+1'; m++)
			C = (quadcross(betaB[|.\m|] , C), 0) + (0, quadcross(kB[|.\m|] , C)) \ C, J(m,1,0)  // prepend next (recursively defined) c polynomial

		C2[zmInd] = colshape(C, 1); C = colshape(C2, `=2*(2*$UAsymParam_M+2)-1') // by shifting right, multiply each c_m(z) poly by z^m; also flip so low-order ones at top

		for (m=2; m<=`=$UAsymParam_M+1'; m++) {
			t1 = *pPT[m]
			t2 = exp(-lngamma(1-beta :: m+1.5-beta))  // exp(lngamma()) faster than gamma()!; even last entry is >0 if nu<0 so lngamma() will return right answer without converting argument to complex
			t3 = t2 :* C[|m,m\m+m,`=3*$UAsymParam_M+1'+m|]
			lnp[m,] = ln(C( (t1 / t2[m  ]) * t3[|.,.\m  ,.|] ))
			lnq[m,] = ln(C( (t1 / t2[m+1]) * t3[|2,.\m+1,.|] ))
		}
		betaDirty = 0
		zDirty = 1
	}

	if (zDirty) {  // evaluate all the p(z) and q(z)
		mlnz =  lnz * (0 .. `=3*$UAsymParam_M+1')
		lnSp = asdfLogRowSumExp(lnp :+ mlnz)
		lnSq = asdfLogRowSumExp(lnq :+ mlnz)
		zDirty = 0
	}

	lnsqrtzdivalpha = 0.5 * (lnz :- lnalpha)
	term = sqrt(*palpha * z); term = C(term + term); lnterm = ln(term)
	return ((1.62e42fefa39efX-001 /*ln(2)*/ + z*.5) :+ ((1 - beta) * lnsqrtzdivalpha - lngammaalpha + asdfLogSumExp(SBessel.lnK(beta-1, term, lnterm) :+ asdfLogSumExp(lnSp :- mlnalpha) \ 
	                                                                                                                SBessel.lnK(beta  , term, lnterm) :+ asdfLogSumExp(lnSq :- mlnalpha) :+ lnsqrtzdivalpha)))
}


class clsStickyFeller {
	real colvector lnBesselKPower_0toM, lnBesselKPower_2Mto0
	complex rowvector lnStickyFeller_tlambdalnu, lnStickyFeller_u2

	real colvector lnStickyFeller()
	complex rowvector lnBesselKPower()  // for consistency, should be given its own class
	void new()
}

void clsStickyFeller::new() {
  complex rowvector u
	u = C(1, (.5..27.5) * 1.6c782f6914edbX-003 /*h*/)  // really (1 + ui)
	lnStickyFeller_u2 = u :* u
	lnStickyFeller_tlambdalnu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ * lnStickyFeller_u2 + ln(u) :- 1.16c0cdab51479X+001 /* + ln 2*2*h/2pi  tack this on here for efficiency */

  lnBesselKPower_0toM = 0::$BesselKPower_M  // stuff for power series
  lnBesselKPower_2Mto0 = lnBesselKPower_0toM * 2
}


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007), Weideman (2010), parabolic contour method to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real colvector t, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu) {
  real colvector _mu, lnp0, negbt, lnattilde, lnm, lnx, lny, nulnx, nulny, zerofill, zeroind, missingfill, lnJacobian, bax, bay; real scalar ba; class clsBesselKAsym SBesselKAsym

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

	real colvector s, xhi, yhi
	real scalar muj, _lnC, lnbax, lnbay, i, j, lngammanegnu, _nulnx, _nulny, beta, _bax, _bay, qa, qb, qc, qd, smaxx, smaxy, Ream1, ta, ta2, ta3, absa2, Ima2, _ix, _iy, absb, anypowerx, anyasymargx, anyasymparamx, anypowery, anyasymargy, anyasymparamy, powerx, asymargx, asymparamx, powery, asymargy, asymparamy
	complex matrix alpha
	complex scalar lnC, lnM, _alpha, lngammaalpha1beta, retval
	complex rowvector lnUAsymParamx, lnUAsymParamy, lnUPowerx, lnUPowery, alphaj
	complex colvector terms, numera
  class clsUPower scalar SUPowerx, SUPowery
  class clsUAsymParam scalar SUAsymParamx, SUAsymParamy

	beta = 1 + nu

	if ((absb=abs(b)) > 1e-10) {
		alpha = _mu * (lnStickyFeller_u2 / -b); if (b > 0) alpha = beta :- alpha  // faster to just compute alpha row by row? (but will still need vessel for retval)

		xhi = bax :> 14.18699
		yhi = bay :> 14.18699
		
		lngammanegnu = lngamma(-nu)
		_lnC = lngamma(beta) - ln(-nu) - lngammanegnu - nu * ln(abs(b)/a)

    SUPowerx = clsUPower(); SUPowerx.setbeta(beta)
    SUPowery = clsUPower(); SUPowery.setbeta(beta)
    SUAsymParamx = clsUAsymParam(); SUAsymParamx.setbeta(beta)
    SUAsymParamy = clsUAsymParam(); SUAsymParamy.setbeta(beta)

		for (j=rows(x);j;j--) {
			lnbax = ln(_bax = bax[j]); _nulnx = nulnx[j]
			lnbay = ln(_bay = bay[j]); _nulny = nulny[j]
      alphaj = alpha[j,]

			// indexes in lambda[j,] after which switch from power or asymptotic-argument evaluation to asymptotic-parameter; min value is 0
			// Spuriously precise thresholds from benchmarking against mpmath's hyperu()
			_iy =  absb / (muj = _mu[j])
      _ix = editmissing(round(sqrt( (xhi[j]? _bax * 2.1846063 : (_bax * .01593272) ^ -1.1629507) * _iy - 1) / 1.6c782f6914edbX-003 /*h*/), 0)
			_iy = editmissing(round(sqrt( (yhi[j]? _bay * 2.1846063 : (_bay * .01593272) ^ -1.1629507) * _iy - 1) / 1.6c782f6914edbX-003 /*h*/), 0)
      anypowerx = xhi[j]==0 & _ix; anyasymargx = xhi[j] & _ix; anyasymparamx = _ix < $StickyFeller_N
      anypowery = yhi[j]==0 & _iy; anyasymargy = yhi[j] & _iy; anyasymparamy = _iy < $StickyFeller_N

			if (anyasymparamx) {
      	SUAsymParamx.setz(_bax, lnbax)
      	SUAsymParamx.setalpha(alphaj[|_ix+1\.|])
        lnUAsymParamx = SUAsymParamx.lnU()
			}
			if (anyasymparamy) {
      	SUAsymParamy.setz(_bay, lnbay)
      	SUAsymParamy.setalpha(alphaj[|_iy+1\.|])
        lnUAsymParamy = SUAsymParamy.lnU()
			}
			if (anypowerx) {
      	SUPowerx.setz(_bax, lnbax)
      	SUPowerx.setalpha(alphaj)
        lnUPowerx = SUPowerx.lnU(_ix)
      }
			if (anypowery) {
      	SUPowery.setz(_bay, lnbay)
      	SUPowery.setalpha(alphaj)
        lnUPowery = SUPowery.lnU(_iy)
      }

      for (i=$StickyFeller_N;i;i--) {
        asymparamx = i > _ix; powerx = anypowerx & asymparamx==0; asymargx = anyasymargx & asymparamx==0
        asymparamy = i > _iy; powery = anypowery & asymparamy==0; asymargy = anyasymargy & asymparamy==0 
				_alpha = alphaj[i]; lngammaalpha1beta = lngamma(_alpha + 1 - beta)
				lnC = _lnC + (lngammaalpha1beta - lngamma(_alpha))

				if (!(powerx & powery))
					lnM = lngammaalpha1beta - lngammanegnu

				if (asymargx | asymargy) {
					Ream1 = Re(_alpha) - 1; absa2 = abs(_alpha); absa2 = absa2*absa2; Ima2 = Im(_alpha); Ima2 = Ima2*Ima2
					ta = Ream1 - nu; ta2 = ta*ta + Ima2; ta3 = -Ream1 - Ream1 - 1 + absa2
					qa = 2*(Ream1+ta); qb = ta3+ta2+4*Ream1*ta; qc = 2*(Ream1*ta2+ta*ta3); qd = ta3*ta2  // coefficients of quartic equation governing when norm of ratio of successive terms is 1
					smaxx = asymargx? max((0, asdfQuarticRoot(qa, qb-_bax*_bax, qc, qd))) : -1 // index about where terms permanently start to increase
					smaxy = asymargy? max((0, asdfQuarticRoot(qa, qb-_bay*_bay, qc, qd))) : -1
					s = 0 :: max((smaxx, smaxy))
					numera = ln(s :+ (_alpha - 1)) + ln(s :+ (_alpha - beta)) - ln(s)  // numerator in large-argument asymptotic series; actually includes part for s! denominator too
				}

				retval = -ln(muj * lnStickyFeller_u2[i] / mu - nu * exp(lnC))
				if (_nulnx < .)
					if (powerx)
						retval = retval + lnUPowerx[i]
					else if (asymparamx)
						retval = retval + lnUAsymParamx[i - _ix] + lnM
					else {  // asymptotic large-argument series
						(terms = (smaxx < smaxy? numera[|.\smaxx+1|] : numera) :- (1.921fb54442d18X+001i /*pi i*/ + lnbax)) [1] =  - _alpha * lnbax
						retval = retval + asdfLogSumExp(quadrunningsum(terms)) + lnM
					}
				if (_nulny < .)
					if (powery)
						retval = retval + lnUPowery[i]
					else if (asymparamy)
						retval = retval + lnUAsymParamy[i - _iy] + lnM
					else {  // asymptotic large-argument series
						(terms = (smaxy < smaxx? numera[|.\smaxy+1|] : numera) :- (1.921fb54442d18X+001i /*pi i*/ + lnbay)) [1] =  - _alpha * lnbay
						retval = retval + asdfLogSumExp(quadrunningsum(terms)) + lnM
					}
				alpha[j,i] = retval  // dicie: sticking return values in alpha even as lnU() routines are holding copies of rows of it for potential reuse
			}
		}

	} else {  // squared Bessel/b=0 case. Would need a different approach if b could vary by observation

		complex rowvector termx, lntermx, termy, lntermy, _lambda, _term, _lnterm, zerolimit; real scalar pidivsinnupi, gammanegnu2, c; real colvector ix, iy

		SBesselKAsym = clsBesselKAsym()
		pidivsinnupi = pi() / sin(nu*pi())
		gammanegnu2 = lngamma(-nu); gammanegnu2 = exp(gammanegnu2 + gammanegnu2)
		zerolimit = J(1, $StickyFeller_N, lngamma(-nu) - ln(2) * beta)
		c =                                            2 * ln(2) * beta 

		iy = (20.25 * a) :/ _mu  // 20.25 = 9^2/4 where 9 is the abs(z) cut-off between power series & asymptotic in BesselK(nu,z)
		ix = editmissing(round( sqrt(iy :/ x :- 1) / 1.6c782f6914edbX-003 /*h*/), 0)  // index of parabola points where BesselK argument passes 9
		iy = editmissing(round( sqrt(iy :/ y :- 1) / 1.6c782f6914edbX-003 /*h*/), 0)

		alpha = J(rows(x), $StickyFeller_N, C(.))
		for (j=rows(alpha);j;j--) {
			_term = sqrt(_lambda = _mu[j] * lnStickyFeller_u2) * (2 / sqrt(a)); _lnterm = ln(_term)  // precompute more of this
			termx = _term * sqrt(x[j]); lntermx = _lnterm :+ .5*lnx[j]
			termy = _term * sqrt(y[j]); lntermy = _lnterm :+ .5*lny[j]

			_ix = ix[j]; _iy = iy[j]
			alpha[j,] = c :+ ((x[j]? (_ix <  1                ? SBesselKAsym.lnK (nu, termx, lntermx) : 
			                         (_ix >= $StickyFeller_N?      lnBesselKPower(nu, lntermx)        :
			                                                       lnBesselKPower(nu, lntermx[|.\_ix|]), SBesselKAsym.lnK(nu, termx[|_ix+1\.|], lntermx[|_ix+1\.|]))) - nu * lntermx :
																zerolimit) +
												(y[j]? (_iy < 1                 ? SBesselKAsym.lnK (nu, termy, lntermy) : 
			                         (_iy >= $StickyFeller_N?      lnBesselKPower(nu, lntermy)        : 
			                                                       lnBesselKPower(nu, lntermy[|.\_iy|]), SBesselKAsym.lnK(nu, termy[|_iy+1\.|], lntermy[|_iy+1\.|]))) - nu * lntermy :
																zerolimit) - 
													 ln(_lambda * (gammanegnu2 / mu) - (_lambda / a) :^ -nu * pidivsinnupi))
		}
	}
	return (asdfLogSumExpRow((lnp0 , lnm + ln(_mu :* rowsum(asdfReexp(lnStickyFeller_tlambdalnu :+ alpha))))))
}


// Compute log Bessel K with power series formula and complex argument z
// series length from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
complex rowvector clsStickyFeller::lnBesselKPower(real scalar nu, complex rowvector lnz) {
  real colvector InuDenom, InegnuDenom; complex matrix twomlnhalfz; complex rowvector nulnhalfz, lnhalfz

  (InuDenom    = ln(( nu::$BesselKPower_M.5+nu) :* lnBesselKPower_0toM)) [1] = lngamma(1+nu)
  (InegnuDenom = ln((-nu::$BesselKPower_M.5-nu) :* lnBesselKPower_0toM)) [1] = lngamma(1-nu)

	lnhalfz = lnz :- 1.62e42fefa39efX-001 /*ln .5*/
	nulnhalfz = nu * lnhalfz
	twomlnhalfz = lnBesselKPower_2Mto0 * lnhalfz
	return ((1.ce6bb25aa1315X-002 /*ln pi/2*/ - ln(C(sin(nu*pi())))) :+ asdfLogSumExp(asdfLogSumExp((twomlnhalfz :- nulnhalfz) :- quadrunningsum(InegnuDenom)) \ 
	                                                                                  asdfLogSumExp((twomlnhalfz :+ nulnhalfz) :- quadrunningsum(   InuDenom)) :+ 1.921fb54442d18X+001i /*pi i*/))
}

complex rowvector asdfLogSumExp(complex matrix x) {
	real rowvector shift
	if (rows(x)==0) return(J(1,cols(x),0))
	if (rows(x)==1) return(x)
	shift = ln(maxdouble()/rows(x)) :- colmax(Re(x))
//	shift = shift - (shift:>0):*shift  // only downshift, to prevent overflow; shifting can prevent underflow & overflow but can also reduce precision if the shifter is much larger than entries
	return (ln(colsum(exp(x :+ shift))) - shift)
}

complex colvector asdfLogRowSumExp(complex matrix x) {
	real colvector shift
	if (cols(x)==0) return(J(rows(x),1,0))
	if (cols(x)==1) return(x)
	shift = ln(maxdouble()/cols(x)) :- rowmax(Re(x))
//	shift = shift - (shift:>0):*shift  // only downshift, to prevent overflow; shifting can prevent underflow & overflow but can also reduce precision if the shifter is much larger than entries
	return (ln(rowsum(exp(x :+ shift))) - shift)
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

* Goran test
mata S = clsStickyFeller(); p = exp(S.lnStickyFeller(t=rangen(.1,2,20)#J(20,1,1), x=J(400,1,1), y=J(20,1,1)#rangen(.1,2,20), a=1, b=1, nu=-.25, mu=5))
drop _all
getmata t y p, force double replace
twoway contour p y t, levels(500) clegend(off)
mata p :/ (exp(b/a*y) :* y:^nu :/ a)
