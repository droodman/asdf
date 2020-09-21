global M 200  // number of terms in power series expansion of M() and U() functions

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

class clsStickyFeller {
	real colvector m, lnm, B

	real colvector lnStickyFeller()
	pointer(real colvector) colvector PT
	complex matrix lnStickyGreenTerm(), lnStickyGreen()
	complex scalar lnUAsymParam(), hypergeomdiff()
	void new()
}

void clsStickyFeller::new() {
	m = 1::$M; lnm = ln(m)

	// B_i / i!, i=1 to 20, B-i = Bernoulli numbers
	B = -1.0000000000000X-001 \ 1.5555555555555X-004 \ 0 \ -1.6c16c16c16c17X-00a \ 0 \ 1.1566abc011566X-00f \ 0 \ -1.bbd779334ef0bX-015 \ 0 \ 1.66a8f2bf70ebeX-01a \ 0 \ -1.22805d644267fX-01f \ 0 \ 1.d6db2c4e09163X-025 \ 0 \ -1.7da4e1f79955cX-02a \ 0 \ 1.355871d652e9eX-02f \ 0 \ -1.f57d968caacf1X-035 \ 0 \ 1.967e1f09c376fX-03a // (-1/2\ 1/6\ 0 \-1/30\ 0 \1/42\ 0 \-1/30\ 0 \5/66\ 0 \-691/2730\ 0 \7/6\ 0 \-3617/510\ 0 \43867/798\ 0 \-174611/330 \ 0 \ 854513/138) :/ factorial(1::22)

	// Pascal's triangle
	PT = &(1\1) \ &(1\2\1) \ &(1\3\3\1) \ &(1\4\6\4\1) \ &(1\5\10\10\5\1) \ &(1\6\15\20\15\6\1) \ &(1\7\21\35\35\21\7\1) \ &(1\8\28\56\70\56\28\8\1) \ &(1\9\36\84\126\126\84\36\9\1) \ &(1\10\45\120\210\252\210\120\45\10\1)
}


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007) parabolic contour method just to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real colvector t, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu) {
  real colvector _mu, lnp0, negbt, lnattilde, lnm, lnx, lny, zerofill, zeroind, missingfill, lnJacobian; complex rowvector u; complex matrix lambda

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
	if (mu == 0)
		return (lnp0)
	if (rows(zerofill)) lnp0[zeroind] = missingfill  // if y=0, "absorbing term" drops out

	_editmissing(lnm = (b / a) * y + nu * lny :- ln(a), -ln(mu))  // log speed measure; takes special value of 1/mu when y = 0

	/*log_epsilon = ln(1e-15)  // Logic extracted from Garrappa's mlf routine for Matlab
	log_eps = ln(epsilon(1))
	w = sqrt(log_eps / (log_eps - log_epsilon))  // half-width of integration range, needed to assure given precision
	N = ceil(-w * log_epsilon / (2*pi())) + 1    // half the number of integration points
	h = w :/ (N-.5) */                           // width of bars in Riemann integral; do an even number of integration points since only integrating over half of parabola

  _mu = 1.813f9e629e9c0X+000 /*log_epsilon - log_eps*/ :/ t  // different "mu" -- as in Weideman and Trefethen eq 1.9
  u = C(1, (.5..27.5) * 1.6c782f6914edbX-003 /*h*/)  // really 1 + ui
  lambda = _mu * (u:*u)
  return (asdfLogSumExpRow((lnp0 , lnm + ln(_mu :* quadrowsum(Re((u * 1.d00e940a4f90fX-004 /* 2*2*h/2pi*/) :* exp(t :* lambda :+ lnStickyGreenTerm(lambda, x, y, a, b, nu, mu, lnx, lny))))))))
}


// log Tricomi U for large a, dlmf.nist.gov/13.8.E11
// optional argument lnM is a multiplier in logs, to be incorporated in order to conserve precision
complex scalar clsStickyFeller::lnUAsymParam(complex scalar a, real scalar b, real scalar z, | complex scalar lnM) {
  real colvector c; real scalar k, s, p, q, t1, t2, t3; complex scalar retval, term, C1, C2, as

	k = 0
	as = p = 1
	c = (q = z/12 - b*.5) \ 1  // c_1 \ c_0 -- store c values backward

	term = sqrt(a * z); term = C(term + term)
	retval = (C1 = BesselK(b-1, term)) + (C2 = q * sqrt(z / a) * BesselK(b, term))
	if (C1==0 & C2==0) return(C(.))
	for (s=1; s<=10; s++) {
		k = k + 2
		t1 = b * B[|1\k+1|] + z * (1::k+1) :* B[|2\k+2|]
		c = c ' t1[|.\k|] /  -k      \ c  // add two more terms to c
		c = c ' t1        / -(k + 1) \ c

		t1 = *PT[s] :* z :^ (s::0)  // since c is stored backward, run these sums s to 0
		t2 = lngamma(s+2-b :: .5-b)
		t3 = c[|1\s+2|] :/ exp(t2)

		p = exp(t2[2]) * quadcross(t1, t3[|2\s+2|])
		q = exp(t2[1]) * quadcross(t1, t3[|1\s+1|])

		retval = retval + (term = (C1 * p + C2 * q) / (as = as * a))
		if (abs(term/retval) < 1e-15) break
	}
	return (1.62e42fefa39efX-001 /*ln(2)*/ + (1 - b)*.5 * ln(z/a) + z*.5 + (editmissing(lnM,0) - lngamma(a)) + ln(retval))
}


// Compute sticky term of sticky Feller Green function for b<0. rows(lambda) must be max(rows(x),rows(y)) and rows(x) and rows(y) must each be that or 1
// accepts optional pre-computed ln(x), ln(y)
complex matrix clsStickyFeller::lnStickyGreenTerm(complex matrix lambda, real colvector x, real colvector y, real scalar a, real scalar b, real scalar nu, real scalar mu, | real colvector lnx, real colvector lny) {
	real colvector t1, t2, denomx1, denomx2, denomy1, denomy2, numer1, numer2, powerxj, poweryj
	complex matrix alpha
	real scalar _lnC, lnbax, lnbay, i, j, N, rowsx, rowsy, lngammanegnu, ba, bax, bay, powers, lnxnu, lnynu, onepnu
	real matrix powerx, powery
	complex scalar lnC, lnM, _alpha, _alphamnu, lngammaalphamnu

	powerxj = rowsum(powerx = abs(lambda :* (x/a)) :< 50)  // which entries to compute using with power series vs. asymptotic
	poweryj = rowsum(powery = abs(lambda :* (y/a)) :< 50)
	alpha = lambda / -b
	lngammanegnu = lngamma(-nu)
	_lnC = lngamma(onepnu = 1 + nu) - ln(-nu) - lngammanegnu - nu * ln(ba = -b/a)

	N = cols(lambda)
	rowsx = rows(x); rowsy = rows(y)
	if (rows(lnx) == 0) lnx = ln(x)
	if (rows(lny) == 0) lny = ln(y)

	t1 = ln(onepnu :: $M.5 + nu) + lnm  // faster than ln(nu :+ m)
	t2 = ln(1 - nu :: $M.5 - nu) + lnm

	if (rowsx == 1) {
		lnxnu = nu * lnx
		bax = ba * x
		if (colsum(powerxj)) {
			lnbax = ln(bax)
			denomx1 = lnbax :- t1
			denomx2 = lnbax :- t2
		}
	}
	
	if (rowsy == 1) {
		lnynu = nu * lny
		bay = ba * y
		if (colsum(poweryj)) {
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

		for (i=N;i;i--) {
			_alpha = alpha[j,i]; _alphamnu = _alpha - nu; lngammaalphamnu = lngamma(_alphamnu)
			lnC = _lnC + (lngammaalphamnu - lngamma(_alpha))
			if (powers = powerx[j,i] + powery[j,i]) {
				numer1 = ln((_alpha    - 1) :+ m)
				numer2 = ln((_alphamnu - 1) :+ m)
			}
			if (powers < 2)
				lnM = lngammaalphamnu - lngammanegnu
			alpha[j,i] = (lnxnu == .? 0 : (powerx[j,i]? hypergeomdiff(numer1 + denomx1, lnC - lnxnu, numer2 + denomx2) : lnUAsymParam(_alpha, onepnu, bax, lnM))) + 
			             (lnynu == .? 0 : (powery[j,i]? hypergeomdiff(numer1 + denomy1, lnC - lnynu, numer2 + denomy2) : lnUAsymParam(_alpha, onepnu, bay, lnM))) -
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
	shift = ln(maxdouble()/rows(x1)) :- max((max(Re(x1)), max(Re(x2))))
	return (ln(exp(shift) - exp(c2 + shift) + quadcolsum(exp(x1 :+ shift) - exp(x2 :+ shift))) - shift)
}


// Compute Bessel K with power series or asymptotic formula for small order (here, -1<nu<1 ASSUMED) and complex argument z
// threshold and index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// except that power series stops at 22 instead of 50, as in worst case (z=9,nu=-.001), this still assures last term is at machine tolerance
// for each row of z, makes a temporary matrix with 23 rows and a column for each row entry with |z|<9, so inefficient if z is very wide
complex matrix BesselK(real scalar nu, complex matrix z) {
  real scalar twonu, K, i, r
	real rowvector k, twokm1, Re_lnz, power, asymp
  real colvector fortyfour2zero, InuDenom, InegnuDenom, lna, z22
  real matrix _power
  complex scalar lnhalfz, nulnhalfz, twomlnhalfz
  complex rowvector _lnz, _z
  complex matrix lnz, terms

  lnz = ln(z)  // *** also using lnz for return value, to avoid allocating twice
  _power = Re(lnz) :< 2.197224577 /* ln 9 */  // abs(z) < 9?

  z22 = 0::22  // stuff for power series
  (InuDenom    = ln(( nu::22.5+nu) :* z22))[1] = lngamma(1+nu); _quadrunningsum(   InuDenom,    InuDenom)
  (InegnuDenom = ln((-nu::22.5-nu) :* z22))[1] = lngamma(1-nu); _quadrunningsum(InegnuDenom, InegnuDenom)
  fortyfour2zero = z22 * 2

  if (nu != -.5 & nu!= -.5) {  // stuff for asymptotic series
  	k = 0::14; twokm1 = 2*k:-1
    twonu = nu + nu
    lna = ln((twonu*twonu :- twokm1:*twokm1) :/ k) :- 1.0a2b23f3bab73X+001 /*ln 8*/  // log of a_k coefficients
  }

  for (r=rows(z);r;r--) {
    power = _power[r,]; asymp = selectindex(1 :- power); power = selectindex(power)
    if (length(power)) {
      lnhalfz = lnz[r,power] :- 1.62e42fefa39efX-001 /*ln .5*/
      nulnhalfz = nu * lnhalfz
      twomlnhalfz = fortyfour2zero * lnhalfz
      lnz[r,power] = 1.921fb54442d18X+000 /*pi/2*/ / sin(nu*pi()) * quadcolsum(exp(twomlnhalfz :- nulnhalfz :- InegnuDenom) - exp(nulnhalfz :+ twomlnhalfz :- InuDenom))
    }
    if (length(asymp))
      if (nu == -.5 | nu == .5) {  // only first term is non-zero
        _z = z[r,asymp]
        lnz[r,asymp] = sqrt(1.921fb54442d18X+000 /*pi/2*/ :/ _z) :/ exp(_z)
      } else {
        _lnz = lnz[r,asymp]
        (terms = lna :- J(15, 1, _lnz))[1,] = 1.ce6bb25aa1315X-003 /*.5*ln(pi/2)*/ :- (.5*_lnz + z[r,asymp])  // log of first term, including outer multiplier, which is sqrt(pi/2z)*exp(-z)
        Re_lnz = Re(_lnz)
        for (i=cols(asymp);i;i--) {
          K = Re_lnz[i]; K = K < 1.c715a530ff3c5X+001 /*ln 35*/ ? 14 : (K < 1.f4bd2b7ac1bafX+001 /*ln 50*/ ? 10 : 8)
          lnz[r,asymp[i]] = quadcolsum(exp(quadrunningsum(terms[|.,i\K,i|], 1)))
        }
      }
  }
  return(lnz)
}

mata mlib create lasdfStickyFeller, dir("`c(sysdir_plus)'l") replace
mata mlib add lasdfStickyFeller *(), dir("`c(sysdir_plus)'l")
mata mlib index
end

mata S = clsStickyFeller()
mata lambda = 149.296316 - 26.7813598i,139.764101 - 80.3440794i,120.699671 - 133.906799i,92.1030264 - 187.469519i,53.9741666 - 241.032238i,6.31309188 - 294.594958i,50.8801978 - 348.157677i,117.605702 - 401.720397i,193.863422 - 455.283117i,279.653356 - 508.845836i,374.975506 - 562.408556i,-479.82987 - 615.971275i,-594.21645 - 669.533995i,718.135244 - 723.096715i,851.586253 - 776.659434i,994.569477 - 830.222154i,1147.08492 - 883.784873i,1309.13257 - 937.347593i,1480.71244 - 990.910313i,1661.82452 - 1044.47303i,1852.46882 - 1098.03575i,2052.64534 - 1151.59847i,2262.35407 - 1205.16119i,2481.59501 - 1258.72391i,2710.36817 - 1312.28663i,2948.67354 - 1365.84935i,3196.51113 - 1419.41207i,3453.88093 - 1472.97479i
mata S.lnStickyGreenTerm(lambda   ,1.5e353f7ced917X-003,1.b4395810624ddX-003,1.1392176df3ec5X+000,-1.42cddd6e04c06X+003,-1.ab7c63fead19bX-002,1.b06a897635e74X+004,-1.c41e964dc9213X+000,-1.8be501affa053X+000)[1]
mata S.lnStickyGreenTerm(lambda[1],1.5e353f7ced917X-003,1.b4395810624ddX-003,1.1392176df3ec5X+000,-1.42cddd6e04c06X+003,-1.ab7c63fead19bX-002,1.b06a897635e74X+004,-1.c41e964dc9213X+000,-1.8be501affa053X+000)

