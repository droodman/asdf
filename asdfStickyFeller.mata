global BesselKPower_M 22  // max term index in power series expansion of BesselK()
global BesselKAsym_M 7    // max term index in asymptotic approximation of BesselK()
global UPower_M 60        // max term index in power series expansion of U()
global UAbadSesma_N  20   // number of terms in Abad-Sesma large-parameter approximation of U()

scalar log_epsilon = ln(1e-15)  // Logic governing evaluation points on parabola of integration, extracted from Garrappa's mlf routine for Matlab
mata st_numscalar("log_eps", ln(epsilon(1)))
scalar _mu = log_epsilon - log_eps  // mu as in Trefethen and Weideman, scaling parabola, not stickiness parameter mu
scalar w = sqrt(-log_eps / _mu)    // half-width of integration range, needed to assure given precision
scalar N = ceil(-w * log_epsilon / (2*_pi)) + 1       //  number of integration points (on positive arm of parabola)
global StickyFeller_h = strofreal(w / (N-.5), "%21x") // width of bars in Riemann integral; do an even number to avoid double-counting evaluation at vertex
if substr("$StickyFeller_h",1,1)=="+" global StickyFeller_h = substr("$StickyFeller_h",2,.)
global StickyFeller_mu = strofreal(_mu, "%21x") // width of bars in Riemann integral; do an even number to avoid double-counting evaluation at vertex
if substr("$StickyFeller_mu",1,1)=="+" global StickyFeller_mu = substr("$StickyFeller_mu",2,.)
global StickyFeller_N = scalar(N)


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

// Compute log (K_nu(z)/z^nu) with asymptotic formula for small order (here, -1<nu<1 ASSUMED) and complex argument z
// max series length index bounds from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
// This code only takes first 8 terms of series, which they recommend for abs(z) > 50.
complex rowvector clsBesselKAsym::lnK(real scalar nu, complex rowvector z, complex rowvector lnz) {
  real scalar twonu, i; complex matrix terms; complex rowvector retval

  if (nu == -.5 | nu == .5)   // only first term is non-zero
		return (1.ce6bb25aa1312X-003 /*ln(sqrt(pi/2))*/ :- ((nu + .5) * lnz + z))

	twonu = nu + nu
	(terms = J(1, i=cols(z), ln(C(twonu*twonu :- twokm1sq)) - k) :- lnz)[1,] = (-.5 - nu) * lnz - z  // log of first term, including outer multipliers
	retval = J(1,i,C(.))
	for (;i;i--)
		retval[i] = 1.ce6bb25aa1312X-003 + asdfLogSumExp(quadrunningsum(terms[,i]))
	return (retval)
}

class clsUPower {
  pointer(complex rowvector) scalar palpha
	real scalar nu, paramDirty, lnz, i
  real colvector m, lnm, mlnz, denom1, denom2, mmnu
  complex colvector terms1, terms2
  complex matrix coefs1, coefs2
  
  void new(), setalpha(), setnu(), setz()
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

void clsUPower::setnu(real scalar _nu) {
  if (nu != _nu) {
    nu = _nu
    paramDirty = 1
    denom1 = quadrunningsum(ln(        nu :: $UPower_M.5 + nu) + lnm)  // first entry will be 0 since lnm[1]=.
    denom2 = quadrunningsum(ln(mmnu = -nu :: $UPower_M.5 - nu) + lnm)
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
complex rowvector clsUPower::lnU(real scalar maxi, complex rowvector lnC) {
  complex scalar _alpha; real scalar _maxi

  _maxi = maxi > $StickyFeller_N? $StickyFeller_N : maxi
	if (paramDirty) i = 1
  for (; i<=_maxi; i++) {
    _alpha = (*palpha)[i]
    (terms1 = ln((_alpha - 1     ) :+ m)) [1] = 0       // for making log Pochhammer symbols (alpha)_m and (alpha+1-beta)_m; edit m=0 entry to ln 1
    (terms2 = ln((_alpha - 1 - nu) :+ m)) [1] = lnC[i]  // edit m=0 entry to multiplier on second series
/*		maxrecoefs[i] = max((max(Re( coefs1[,i] = quadrunningsum(terms1) - denom1 )),
		                     max(Re( coefs2[,i] = quadrunningsum(terms2) - denom2 ))))
"_alpha"
_alpha
"quadrunningsum(terms1), (lngamma((_alpha -  1      ) :+ m)) :-  lngamma(_alpha - 1     )"
 quadrunningsum(terms1), (lngamma((_alpha -  1      ) :+ m)) :-  lngamma(_alpha - 1     )*/
		coefs1[,i] =  (lngamma( _alpha       :+ m) :- lngamma(_alpha     )) - denom1
		coefs2[,i] = ((lngamma((_alpha - nu) :+ m) :- lngamma(_alpha - nu)) - denom2) :+ lnC[i]
  }
  paramDirty = 0

/*	shift = (lnz > 0? `=ln(1.fffffffffffffX+3fe/2/($UPower_M+1))' - ($UPower_M - nu) * lnz : `=ln(1.fffffffffffffX+3fe/2/($UPower_M+1))') :- maxrecoefs
	return (editmissing(ln(quadcolsum(exp(coefs1 :+       mlnz :+ shift)) - 
	                       quadcolsum(exp(coefs2 :+ mmnu * lnz :+ shift))) - shift, 0)) */
/*`"strofreal(Re(exp(coefs1 :+ mlnz)),"%30.15f") :+ " + " :+ strofreal(Im(exp(coefs1 :+ mlnz)),"%30.15f") :+"j""'
  strofreal(Re(exp(coefs1 :+ mlnz)[,1]),"%30.15g") :+ " + " :+ strofreal(Im(exp(coefs1 :+ mlnz)[,1]),"%30.15g") :+"j",strofreal(Re(quadrunningsum(exp(coefs1 :+ mlnz)[,1])),"%30.15g") :+ " + " :+ strofreal(Im(quadrunningsum(exp(coefs1 :+ mlnz)[,1])),"%30.15g") :+"j"*/
/*"(coefs1 :+ mlnz \ coefs2 :+ (1.921fb54442d18X+001i :+ mmnu * lnz))[,6],ln(quadrunningsum(exp(coefs1 :+ mlnz \ coefs2 :+ (1.921fb54442d18X+001i :+ mmnu * lnz))[,6]))"
 (coefs1 :+ mlnz \ coefs2 :+ (1.921fb54442d18X+001i :+ mmnu * lnz))[,6],ln(quadrunningsum(exp(coefs1 :+ mlnz \ coefs2 :+ (1.921fb54442d18X+001i :+ mmnu * lnz))[,6]))*/
 
	return (editmissing(asdfLogSumExp(coefs1 :+ mlnz \ coefs2 :+ (1.921fb54442d18X+001i :+ mmnu * lnz)), 0))  // return power series sum in logs; lnU(0) = 0
}


class clsUAsymArg {
	real scalar alphaDirty, nuDirty, nu, z, lnz, lngammanegnu
	real rowvector qa, qb, qc, qd
	complex rowvector lngammaalphamnu
	pointer(complex rowvector) scalar palpha

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
// fails for entries of a whose real part is a negative integer because then lngamma(a) wrongly returns missing (bug reported to Stata 9/29/20)
complex rowvector clsUAsymArg::lnU(real scalar maxi) {
	real rowvector Ream1, ta, ta2, ta3, absa2, Ima2, M; complex rowvector retval; real colvector m, numer; complex colvector terms; real scalar i, _maxi

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
	for (i=cols(M);i;i--) {
		m = 0 :: M[i]
		numer = ln(m :+ ((*palpha)[i] - 1)) + ln(m :+ ((*palpha)[i] - 1 - nu)) - ln(m)  // scope for saving and reusing these when only z changes
		(terms = numer :- (1.921fb54442d18X+001i /*pi i*/ + lnz)) [1] =  lngammaalphamnu[i] - lngammanegnu - (*palpha)[i] * lnz  // scope for simplification if move from runningsum(ln()) to lngamma()
		retval[i] = asdfLogSumExp(quadrunningsum(terms))
	}
	alphaDirty = nuDirty = 0
	return (retval)
}


// log Tricomi U for large a, using Abad and Sesma (1997). Normalized so U(a,b,0)=1
class clsUAbadSesma {
  real scalar nu, z, lnz, nuDirty, beta, onemhalfbeta
	real colvector B, f
  complex colvector g
  pointer(complex rowvector) scalar palpha
  class clsBesselKAsym scalar S

  void new(), setalpha(), setnu(), setz()
	complex rowvector lnU()
}

void clsUAbadSesma::new() {
// mata B = J(100,1,1)
// python: import mpmath as mp
// forvalues i=1/100 {
//   python:  Mata.storeAt("B", `i'-1, 0, 2**(2*`i')*abs(mp.bernoulli(2*`i')/`i'))
// }
// mata strofreal(B, "%21x")

  // 4^i/i*|B_2i|, B Bernoulli numbers; first element is for i=1
  B = 1.5555555555555X-001\1.1111111111111X-002\1.0410410410410X-001\1.1111111111111X+001\1.f07c1f07c1f08X+003\1.5995995995995X+007\1.5555555555556X+00b\1.c5e5e5e5e5e5eX+00f\1.86e7f9b9fe6e8X+014\1.a74ca514ca515X+019\1.1975cc0ed7304X+01f\1.c2f0566566567X+024\1.ac572aaaaaaabX+02a\1.dc0b1a5cfbe17X+030\1.31fad7cbf3c00X+037\1.c280563b8bcbdX+03d\1.7892edfdf5556X+044\1.62b8b44651d09X+04b\1.76024c215d22bX+052\1.b6c0dfed2955bX+059\1.1cca39b77b027X+061\1.97212d8cc1040X+068\1.3f0cb06b17e29X+070\1.1101d96823ee1X+078\1.fc474bdd53c20X+07f\1.007db56db95dfX+088\1.17c6dd28a9378X+090\1.48df88a383ad8X+098\1.9f7b3fa37f314X+0a0\1.195c16c40d563X+0a9\1.97922eafb5d17X+0b1\1.3b0a43def5904X+0ba\1.035a171273533X+0c3\1.c5e89f1fd242eX+0cb\1.a575dd47b788aX+0d4\1.9e84a01ae153bX+0dd\1.af26959307253X+0e6\1.d9890bc68eeddX+0ef\1.1231831c6ccbdX+0f9\1.4e5b1c79769adX+102\1.acc2917790916X+10b\1.20bd5e935fc86X+115\1.97f9f572293ebX+11e\1.2e083516087a9X+128\1.d41e650690e84X+131\1.7b59a0f58ffb8X+13b\1.412709736939fX+145\1.1bc48bbc5cc50X+14f\1.0575ec3ab2234X+159\1.f5fe45963cac3X+162\1.f5abf07320dbaX+16c\1.04c0734d5322dX+177\1.19bdb7ca080b5X+181\1.3c2f25da2f119X+18b\1.704a6714a172bX+195\1.bcf2002532ee2X+19f\1.169855a5defafX+1aa\1.6963e57ed41dfX+1b4\1.e54fbbb872117X+1be\1.5125bc119ad89X+1c9\1.e4616e3d8d09aX+1d3\1.679a15f64d5a9X+1de\1.13c15fcb28309X+1e9\1.b49d69397c7aaX+1f3\1.64abed86e0175X+1fe\1.2c8197ead1271X+209\1.05015c2200ecfX+214\1.d32f2c47ba10cX+21e\1.aea4c26df15f6X+229\1.98adca0fa6ca5X+234\1.8f1dcfb3cd633X+23f\1.90f42af6fe074X+24a\1.9e2f8a971c6f1X+255\1.b7c7b687b5bb6X+260\1.dfcb00aa4c2cfX+26b\1.0cd20cf8a2e4dX+277\1.354a937fa88e0X+282\1.6d400c144ce9eX+28d\1.ba92b59047eb7X+298\1.13079715a5f66X+2a4\1.5e8194cdc3f82X+2af\1.c9e5d8fd49582X+2ba\1.32811f8ce591dX+2c6\1.a45eae570d41eX+2d1\1.273c29596492dX+2dd\1.a89a664a7f630X+2e8\1.3888d893937f4X+2f4\1.d6d21fe9b4fb0X+2ff\1.6acfec47defa2X+30b\1.1df41024fe6a3X+317\1.cce916a2f02f1X+322\1.7bbbe1e92e94bX+32e\1.3fbfc0ec34e1dX+33a\1.131bca925b254X+346\1.e39bd06647752X+351\1.b221ac47c86e0X+35d\1.8df13eed18d1cX+369\1.746394821e395X+375\1.63ae13ac3b025X+381\1.5aabe04f8c67bX+38d

  f = J($UAbadSesma_N,1,1)
  g = J($UAbadSesma_N,1,C(1))
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

complex rowvector clsUAbadSesma::lnU(real scalar mini) {
  real colvector j; complex scalar Karg, lnKarg, lnnegiz; real scalar s; complex matrix terms

  if (nuDirty) {
  	beta = 1 + nu
    onemhalfbeta = 1 - .5*beta
    for (s=1;s<$UAbadSesma_N;s++) {
      j = 0 :: s-1
      f[s+1] = onemhalfbeta * quadcross(comb(2*s-1, 2*j), B[s:-j], f[|.\s|])
    }
    nuDirty = 0
  }

  lnnegiz = -1.921fb54442d18X+000i /*log -i*/ + lnz

  Karg = sqrt((4 * (*palpha)[|mini+1\.|] :- 2 * beta) * z); lnKarg = ln(Karg)
  (terms = J($UAbadSesma_N, cols(Karg), C(.))) [1,] = S.lnK(beta-1, Karg, lnKarg)
  for (s=1;s<$UAbadSesma_N;s++) {
    j = 0 :: floor((s - 1) * .5)
    g[s+1] = -.25i * (z * ((comb(s-1, 2*j) :* B[j:+1]) ' g[s:-2*j]))
    j = 0 :: floor(s * .5)
    terms[s+1,] = (s * lnnegiz - lngamma(s+1) + ln((comb(s, 2*j) :* f[j:+1]) ' g[(s+1):-2*j])) :+ S.lnK(beta-1+s, Karg, lnKarg)  // small bug: lnK not written for order N+.5 other than +-.5
  }

  return (asdfLogSumExp(terms) :+ (beta * 1.62e42fefa39efX-001 /*ln(2)*/ + z*.5 - lngamma(1 - beta)))
}


class clsStickyFeller {
	real colvector lnBesselKPower_0toM, lnBesselKPower_2Mto0, ot, ox, oy, zerofill, zeroind, missingfill, lnx, lny, Nuniq, _mu
	real matrix uniqdata
	pointer (real colvector) scalar px, py, pt
	real scalar N, Nt
	complex rowvector tlambdalnu, u2, zeros
	complex colvector phi
	complex matrix lambda
  class clsUPower scalar SUPower  // instantiate these objectsion when clsStickyFeller is instantiated--just once per estimation
  class clsUAbadSesma scalar SUAsymParam
  class clsUAsymArg scalar SUAsymArg
	class clsBesselKAsym scalar SBesselKAsym

	real colvector lnStickyFeller()
	complex rowvector lnBesselKPower()  // for consistency, should be given its own class
	void new(), setData()
}

void clsStickyFeller::new() {
  complex rowvector u
	
  lnBesselKPower_0toM = 0::$BesselKPower_M  // stuff for power series
  lnBesselKPower_2Mto0 = lnBesselKPower_0toM * 2
	
	u = C(1, (.5 .. $StickyFeller_N-.5) * $StickyFeller_h)  // really (1 + ui)
	u2 = u :* u
	tlambdalnu = ln(u) + $StickyFeller_mu * u2 :`=strofreal(ln(2*2*$StickyFeller_h/2/_pi), "%21x")' // + ln 2*2*h/2pi for Mellin integral; tack on here for efficiency. strofreal() always includes +-
	
	zeros = J(1, $StickyFeller_N, C(0))
}


void clsStickyFeller::setData(real colvector t, real colvector x, real colvector y) {
	real colvector tID, txyID, o, invo; real matrix data, tinfo, txyinfo
	pragma unset tID; pragma unset txyID

	px = &x; py = &y; pt = &t
	
	zerofill = J(length(zeroind = selectindex(y:==0)), 1, 0)
	missingfill = J(length(zeroind), 1, .)
	lnx = ln(x); lny = ln(y)

	N = rows(t)
	data = t,x,lnx \ t,y,lny
	_collate(data, o = order(data, (1,2)))  // sort observations by t, x/y
	tinfo   = asdfpanelsetup(data,  1   ,   tID)'
	txyinfo = asdfpanelsetup(data, (1,2), txyID)
	uniqdata = data[txyinfo[,1],]  // unique t-x/y combinations
	Nuniq = rows(txyinfo)
	Nt = cols(tinfo)  // unique values of t
	invo = invorder(o)
	ot =   tID[invo[|.  \N|]]  // map from original observations to unique t's
	ox = txyID[invo[|.  \N|]]  // map from original observations to unique t-x's
	oy = txyID[invo[|N+1\.|]]  // map from original observations to unique t-y's
	phi = J(Nuniq, $StickyFeller_N, C(.))

	_mu = $StickyFeller_mu :/ data[tinfo[1,],1]
	lambda = _mu * u2   // "mu" as in Garrappa (2015), p. 1364, not stickiness parameter
}


// compute log transition density f(t; x, y) = p(t; x, y) * m(dx) for sticky Feller by applying Weideman & Trefethen (2007), Weideman (2010), parabolic contour method to sticky term of Green function
// assumes b < 0 and t,x,y have same height
real colvector clsStickyFeller::lnStickyFeller(real scalar a, real scalar lna, real scalar b, real scalar nu, real scalar mu, real scalar lnmu) {
  real colvector lnp0, negbt, lnattilde, lnm, lnJacobian, z, lnz
	real scalar hi, ba, i, j, lngammanegnu, beta, p, absb, anypower, anyasymparam, anyasymarg, oldt, lnba, absbdiv_mu, c, floorp
	real rowvector zerolimit
	complex matrix alpha, lngammaalphamnu, K, lnC
	complex rowvector lnUAsymParam, lnUPower, lnUAsymArg, alphaj, lnalphaj, lambdaj, lnCj
	pragma unset oldt  // will exploit that initial value is missing
  pragma unset zerolimit

	if (nu < -1 | nu > 0)
		return (.)

	if (b)
		lnattilde = ln(a / (-b) * expm1(negbt = (-b) * *pt))
	else {
		lnattilde = ln(a * *pt)
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

	beta = 1 + nu
	lngammanegnu = lngamma(-nu)

	absb = abs(b)
	ba = absb / a
	_editmissing(lnm = nu * lny - ba * (b < 0? *py : *px) :- lna, -lnmu)  // log speed measure; takes special value of 1/mu when y = 0; b>0 modification is a hack to multiply by exp(-b/a*(x+y))

	if (absb > 1e-10) {
		SUPower.setnu    (nu)
		SUAsymParam.setnu(nu)
		SUAsymArg.setnu  (nu)

		alpha = (lambda / -b); if (b > 0) alpha = beta :- alpha

		lngammaalphamnu = lngamma(alpha :- nu)
		lnC = (lngamma(beta) - ln(-nu) - lngammanegnu) :+ (lngammaalphamnu - lngamma(alpha))
		K = -ln(u2 / mu :- (nu :/ _mu) :* exp(lnC :- nu * (lnba = ln(absb)-lna)))  // K in G = K * phi(x) * phi(y); mu is stickiness parameter; _mu is mu as in Garrappa (2015), p. 1364 (here, part of Jacobian of parabolic path function)

		j = Nt + 1
		for (i=Nuniq; i; i--) {
			if (uniqdata[i,1] != oldt) {
				oldt = uniqdata[i,1]
				alphaj = alpha[--j,]
				absbdiv_mu = absb / _mu[j]
				lnCj = .
				SUPower.setalpha(alphaj)
				SUAsymParam.setalpha(alphaj)
				SUAsymArg.setalpha(alphaj)
			}
			z = ba * uniqdata[i,2]; lnz = lnba :+ uniqdata[i,3]

			// indexes in lambda after which switch from power or asymptotic-argument  to asymptotic-parameter; min value is 0
			// Spuriously precise thresholds from benchmarking against mpmath's hyperu()
			hi = z :& z :> 15.460074  // in U() input region where low-z values best approximated by large-argument asymptotic rather than standard power series representation?
			p = editmissing(.5 + sqrt( (hi? (z * `=exp(-1.446009 /  1.102413)') ^  1.102413 : 
                                      (z * `=exp( 3.582745 / -.7823305)') ^ -.7823305   ) * absbdiv_mu - 1) / $StickyFeller_h, 0)  // index of last point on parabola handled with power series representation
			floorp = floor(p)
      anypower = hi==0 & p >= 0; anyasymarg = hi & p >= 0; anyasymparam = p < $StickyFeller_N

			if (anypower) {
				SUPower.setz(z, lnz)
				lnUPower = SUPower.lnU(floorp+1, (lnCj==.? (lnCj=lnC[j,]) : lnCj))
			} else if (anyasymarg) {
				SUAsymArg.setz(z, lnz)
				lnUAsymArg = SUAsymArg.lnU(floorp+1)
			}
			if (anyasymparam) {
				SUAsymParam.setz(z, lnz)
				lnUAsymParam = SUAsymParam.lnU(floorp)
			}
      if (p >= 0 & anyasymparam) {  // if using both power series/large-argument and large-parameter across lambda vector, smooth cross-over with a weighted average to reduce discontinuities w.r.t parameters
        c = p - floorp
        lnUAsymParam[1] = (1 - c) *lnUAsymParam[1] + c * (hi? lnUAsymArg : lnUPower)[floorp+1]
      }
			phi[i,] = z? (hi? (p >= 1? (p >= $StickyFeller_N? lnUAsymArg : lnUAsymArg[|.\floorp|], lnUAsymParam) : lnUAsymParam) :  // value of eigenfunction factor for this t and x/y
			                  (p >= 1? (p >= $StickyFeller_N? lnUPower   : lnUPower  [|.\floorp|], lnUAsymParam) : lnUAsymParam))      :
										zeros
/*if (i==46) {
"i,hi,z,nu,p,anypower,anyasymarg,anyasymparam,"
 i,hi,z,nu,p,anypower,anyasymarg,anyasymparam
"alphaj\phiy\lnC[j,]\SUPower.lnU($StickyFeller_N, lnC[j,])\SUAsymParam.lnU(0)"
 alphaj\phi[i,]\lnC[j,]\SUPower.lnU($StickyFeller_N, lnC[j,])\SUAsymParam.lnU(0)
"z,RelnCj[6],ImlnCj[6],Realpha[6], Imalpha[6]"
strofreal((z,Re(lnCj[6]),Im(lnCj[6])),"%21x"),strofreal((Re(alphaj[6]),Im(alphaj[6])),"%21x")
external numeric matrix _phi, _alpha, _lnC, _z, _nu
_phi=phi[i,];_alpha=alphaj; _lnC=lnCj; _z=z; _nu=nu
if (abs(nu- -.096)<1e-6 & 0)
  I(1)[2]
}*/

		}

	} else {  // squared Bessel/b=0 case. Would need a different approach if b could vary by observation
		lnC = (lngamma(beta) - ln(-nu) - lngammanegnu) :- nu * ln(lambda / a)
		K = -ln(u2 / mu :- (nu :/ _mu) :* exp(lnC)) :-   // K in G = K * phi(x) * phi(y); mu is stickiness parameter; _mu is mu as in Garrappa (2015), p. 1364 (here, part of Jacobian of parabolic path function)
		                                    nu * 1.62e42fefa39efX+000 /*2*ln(2)*/  // offsets multiplying phi(x) and phi(y) by 2^nu to give them computationally efficient form X*K_nu(X)
		j = Nt + 1
		for (i=Nuniq; i; i--) {
			if (uniqdata[i,1] != oldt) {
				oldt = uniqdata[i,1]
				lambdaj = lambda[--j,]
			}

			z = uniqdata[i,2]
			
			alphaj = sqrt(lambdaj :* ((4/a) * z)); lnalphaj = ln(alphaj)

			p = editmissing(round( sqrt((20.25 * a / _mu[j]) / z :- 1) / 1.6c782f6914edbX-003 /*h*/), 0)  // indexes of parabola points where BesselK argument passes 9, cut-off between power series & asymptotic; 20.25 = 9^2/4 

			phi[i,] = z? (p <  1 ?                                                                   SBesselKAsym.lnK(nu, alphaj,          lnalphaj         )  : 
			                        (p >= $StickyFeller_N?     lnBesselKPower(nu, lnalphaj)        :
			                                                   lnBesselKPower(nu, lnalphaj[|.\p|]), SBesselKAsym.lnK(nu, alphaj[|p+1\.|], lnalphaj[|p+1\.|]))) :
									           (cols(zerolimit)? zerolimit : (zerolimit = J(1, $StickyFeller_N, nu * 1.62e42fefa39efX-001 /*ln(2)*/)))  // real z->0 limit is 1, but we're multiplying everything by 2^nu so phi has form X*K_nu(X)
		}
	}
/*external complex matrix _tlambdalnu, _K, _phix, _phiy, _lnC
_tlambdalnu=_tlambdalnu\tlambdalnu; _K=_K\K[ot,][373,]; _phix=_phix\phi[ox,][373,]; _phiy=_phiy\phi[oy,][373,]; _lnC=lnCj
"tlambdalnu\ K[ot,][373,]\ phi[ox,][373,] \ phi[oy,][373,]"
 tlambdalnu\ K[ot,][373,]\ phi[ox,][373,] \ phi[oy,][373,]
"ot[373],ox[373],oy[373],lnp0[373],lnm[373],asdfLogRowSumExpC((tlambdalnu :+ K)[ot,] + phi[ox,] + phi[oy,])[373]"
 ot[373],ox[373],oy[373],lnp0[373],lnm[373],asdfLogRowSumExpC((tlambdalnu :+ K)[ot,] + phi[ox,] + phi[oy,])[373]*/

	return (Re(asdfLogRowSumExpC((lnp0 , lnm + asdfLogRowSumExpC((tlambdalnu :+ K)[ot,] + phi[ox,] + phi[oy,])))))  // In princple can take Re() sooner, but in extreme cases imprecision makes Mellin "integral" <0, throwing us into C
}


// Compute log (K_nu(z) / z^nu) with power series formula and complex argument z
// series length from Zhang and Jin's CIKVB @ people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90
complex rowvector clsStickyFeller::lnBesselKPower(real scalar nu, complex rowvector lnz) {
  real colvector InuDenom, InegnuDenom; complex matrix twomlnhalfz; complex rowvector lnhalfz

  (InuDenom    = ln(( nu::$BesselKPower_M.5+nu) :* lnBesselKPower_0toM)) [1] = lngamma(1+nu)
  (InegnuDenom = ln((-nu::$BesselKPower_M.5-nu) :* lnBesselKPower_0toM)) [1] = lngamma(1-nu)

	lnhalfz = lnz :- 1.62e42fefa39efX-001 /*ln .5*/
	twomlnhalfz = lnBesselKPower_2Mto0 * lnhalfz
	return ((1.ce6bb25aa1315X-002 /*ln pi/2*/ - nu * 1.62e42fefa39efX-001 /*ln 2*/ - ln(C(sin(nu*pi())))) :+ 
                asdfLogSumExp(asdfLogSumExp((twomlnhalfz :+ nu * (1.62e42fefa39efX-000 /*ln 4*/ :- 2*lnz)) :- quadrunningsum(InegnuDenom)) \ 
	                            asdfLogSumExp( twomlnhalfz                                                   :- quadrunningsum(   InuDenom)) :+ 1.921fb54442d18X+001i /*pi i*/))
}
mata
complex rowvector asdfLogSumExp(complex matrix x) {
	real rowvector shift; real colvector limits, minn, maxx; real matrix t
	if (rows(x)==0) return(J(1,cols(x),0))
	if (rows(x)==1) return(x)
	limits = ln(smallestdouble()) + 1 \ ln(maxdouble()) - rows(x) - 1
	t = limits :- colminmax(Re(x))  // shift just enough to prevent underflow in exp(), but if necessary even less to avoid overflow
	minn = t[1,]; maxx = t[2,]
	t = minn :* (minn :> 0) - maxx
	shift = t :* (t :< 0) + maxx // parallelizes better than rowminmax()
	return (any(shift)? ln(quadcolsum(exp(x :+ shift))) - shift : ln(quadcolsum(exp(x))))
}

complex colvector asdfLogRowSumExp(complex matrix x) {
	real colvector shift
	if (cols(x)==0) return(J(rows(x),1,0))
	if (cols(x)==1) return(x)
	shift = ln(maxdouble()/cols(x)) :- rowmax(Re(x))
	return (ln(rowsum(exp(x :+ shift))) - shift)
}

// sum pairs of numbers stored in logs, avoiding overflow, treating missing as log of 0
complex colvector asdfLogRowSumExpC(complex matrix x) {
	real colvector shift; real rowvector limits, minn, maxx; real matrix t
	limits = ln(smallestdouble()) + 1 , ln(maxdouble()) - rows(x) - 1
	t = limits :- rowminmax(Re(x))  // shift just enough to prevent underflow in exp(), but if necessary even less to avoid overflow
	minn = t[,1]; maxx = t[,2]
	t = minn :* (minn :> 0) - maxx
	shift = t :* (t :< 0) + maxx // parallelizes better than rowminmax()
	return (any(shift)? ln(quadrowsum(exp(x :+ shift))) - shift : ln(quadrowsum(exp(x))))
}


// Re(exp(z)) but a bit faster
real matrix asdfReexp(complex matrix z)	return (exp(Re(z)) :* cos(Im(z)))

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

mata mlib create lasdfStickyFeller, dir("`c(sysdir_plus)'l") replace
mata mlib add lasdfStickyFeller *(), dir("`c(sysdir_plus)'l")
mata mlib index
end

cap mata mata drop lnUmpmath()
mata
complex scalar lnUmpmath(alpha, beta, z) {
  external _alpha, _beta, _z, _t
  _alpha=alpha; _beta=beta; _z=z; _t = C(.)
  stata(`"python: alpha = Mata.getAt("_alpha",0,0)"')
  stata(`"python: beta  = Mata.getAt("_beta",0,0)"')
  stata(`"python: z     = Mata.getAt("_z",0,0)"')
  stata(`"cap python: Mata.store("_t", complex(mp.ln(mp.hyperu(alpha, beta, z))-mp.ln(mp.hyperu(alpha, beta, 0))))"')
  return(_t)
}
end

python: import mpmath as mp
python: from sfi import *
mata S=clsUPower();S.setalpha(J(1,28,1.380c75879bc19X+002+1.c70bd809591aaX+007i));S.setz(1.00df729aa5153X-001,ln(1.00df729aa5153X-001));lnC=J(1,28,1.43e07c4625a9cX-001+1.310248f9e3000X-003i);S.setnu(-.097);S.lnU(28,lnC)[6];S.setnu(-.096);S.lnU(28,lnC)[6];S.setnu(-.095);S.lnU(28,lnC)[6]
mata S=clsUAbadSesma();S.setalpha(J(1,28,1.380c75879bc19X+002+1.c70bd809591aaX+007i));S.setz(1.00df729aa5153X-001,ln(1.00df729aa5153X-001));lnC=J(1,28,1.43e07c4625a9cX-001+1.310248f9e3000X-003i);S.setnu(-.097);S.lnU(0)[6];S.setnu(-.096);S.lnU(0)[6];S.setnu(-.095);S.lnU(0)[6]
mata lnUmpmath(1.380c75879bc19X+002+1.c70bd809591aaX+007i, 1-.097, 1.00df729aa5153X-001),lnUmpmath(1.380c75879bc19X+002+1.c70bd809591aaX+007i, 1-.096, 1.00df729aa5153X-001),lnUmpmath(1.380c75879bc19X+002+1.c70bd809591aaX+007i, 1-.095, 1.00df729aa5153X-001)
---
mata
"Not so accurate. mpmath and Wolfram agree on -9.22389593194306062766860063988491893259125722132768440312874552... + 2.33503404129907674957498389809945046830061346887716225232505179... i"
alpha=-310.1341902724561805371+132.2627653666411049471i;beta=0.9552103; z=0.6249443861894480312; nu=beta-1
lnC = lngamma(beta) - ln(-nu) - lngamma(-nu) + lngamma(alpha-nu) - lngamma(alpha)
S=clsUAsymParam();S.setnu(nu);S.setalpha(alpha);S.setz(z,ln(z));S.lnU(0)
S=clsUPower    ();S.setnu(nu);S.setalpha(alpha);S.setz(z,ln(z));S.lnU(1,lnC)[1]
S=clsUAbadSesma();S.setnu(nu);S.setalpha(alpha);S.setz(z,ln(z));S.lnU(0)
end

---
* Goran test
mata S = clsStickyFeller(); S.setData(t=rangen(.1,2,20)#J(20,1,1), x=J(20*20,1,1), y=J(20,1,1)#rangen(.1,2,20)); p = exp(S.lnStickyFeller(a=exp(1, ln(a), b=-1, nu=-.25, mu=5, ln(mu))) // pretty: (a=.2, ln(a), b=-1, nu=-.75, mu=5, ln(mu))
drop _all
getmata t y p, force double replace
twoway contour p y t, levels(500) clegend(off) plotregion(margin(zero)) scheme(s1rcolor)
mata quadcolsum(p :/ (exp(b/a*y) :* y:^nu :/ a))

