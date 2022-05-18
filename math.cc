#include <cstdint>

#define _FOFF  7
#define _FFRAC  ((unsigned short)((1 << _FOFF) - 1))
#define _FMASK  ((unsigned short)(0x7fff & ~_FFRAC))
#define _FMAX   ((unsigned short)((1 << (15 - _FOFF)) - 1))
#define _FSIGN  ((unsigned short)0x8000)

#define _DENORM    (-2)
#define _FINITE    (-1)
#define _INFCODE   1
#define _NANCODE   2

union _Fval
{
	unsigned short _Sh[ 8 ];
	float _Val;
};

#define _F0 1 // little-endian
#define _F1 0

extern "C" short _fdtest( float* px )
{
	const auto ps = reinterpret_cast<_Fval*>( px );

	if ( ( ps->_Sh[ _F0 ] & _FMASK ) == _FMAX << _FOFF )
	{
		return static_cast<short>( ( ps->_Sh[ _F0 ] & _FFRAC ) != 0 || ps->_Sh[ _F1 ] != 0 ? _NANCODE : _INFCODE );
	}
	else if ( ( ps->_Sh[ _F0 ] & ~_FSIGN ) != 0 || ps->_Sh[ _F1 ] != 0 )
	{
		return ( ps->_Sh[ _F0 ] & _FMASK ) == 0 ? _DENORM : _FINITE;
	}
	else
	{
		return 0;
	}
}

extern "C" int _fltused = 0;

union float_long
{
	float f;
	long l;
};

float frexpf( float x, int* pw2 )
{
	union float_long fl;
	long int i;

	fl.f = x;
	/* Find the exponent (power of 2) */
	i = ( fl.l >> 23 ) & 0x000000ff;
	i -= 0x7e;
	*pw2 = i;
	fl.l &= 0x807fffff; /* strip all exponent bits */
	fl.l |= 0x3f000000; /* mantissa between 0.5 and 1 */
	return( fl.f );
}

float ldexpf( float x, int pw2 )
{
	union float_long fl;
	long e;

	fl.f = x;

	e = ( fl.l >> 23 ) & 0x000000ff;
	e += pw2;
	fl.l = ( ( e & 0xff ) << 23 ) | ( fl.l & 0x807fffff );

	return( fl.f );
}

extern "C" float sqrtf( const float x )
{
	float f, y;
	int n;

	if ( x == 0.0 ) return x;
	else if ( x == 1.0 ) return 1.0;
	else if ( x < 0.0 ) return 0.0;

	f = frexpf( x, &n );
	y = 0.41731 + 0.59016 * f; /*Educated guess*/
	/*For a 24 bit mantisa (float), two iterations are sufficient*/
	y += f / y;
	y = ldexpf( y, -2 ) + f / y; /*Faster version of 0.25 * y + f/y*/

	if ( n & 1 )
	{
		y *= 0.7071067812;
		++n;
	}
	return ldexpf( y, n / 2 );
}

static float FOPI = 1.27323954473516;

static float PIF = 3.14159265358979323846F;
static float PIO4F = 0.7853981633974483096;
static float PIO2F = 1.570796326794896619F;

/* These are for a 24-bit significand: */
static float DP1 = 0.78515625;
static float DP2 = 2.4187564849853515625e-4;
static float DP3 = 3.77489497744594108e-8;
static float lossth = 8192.;
static float T24M1 = 16777215.;

static float sincof[] = {
	-1.9515295891E-4,
	8.3321608736E-3,
	-1.6666654611E-1
};
static float coscof[] = {
	2.443315711809948E-005,
	-1.388731625493765E-003,
	4.166664568298827E-002
};

extern "C" float sinf( float xx )
{
	float* p;
	float x, y, z;
	/*register*/ unsigned long j;
	/*register*/ int sign;

	sign = 1;
	x = xx;
	if ( xx < 0 )
	{
		sign = -1;
		x = -xx;
	}
	if ( x > T24M1 )
	{
		// mtherr( "sinf", TLOSS );
		return( 0.0 );
	}
	j = FOPI * x; /* integer part of x/(PI/4) */
	y = j;
	/* map zeros to origin */
	if ( j & 1 )
	{
		j += 1;
		y += 1.0;
	}
	j &= 7; /* octant modulo 360 degrees */
			/* reflect in x axis */
	if ( j > 3 )
	{
		sign = -sign;
		j -= 4;
	}

	if ( x > lossth )
	{
		// mtherr( "sinf", PLOSS );
		x = x - y * PIO4F;
	}
	else
	{
		/* Extended precision modular arithmetic */
		x = ( ( x - y * DP1 ) - y * DP2 ) - y * DP3;
	}
	/*einits();*/
	z = x * x;
	if ( ( j == 1 ) || ( j == 2 ) )
	{
		/* measured relative error in +/- pi/4 is 7.8e-8 */
		/*
		y = ((  2.443315711809948E-005 * z
		- 1.388731625493765E-003) * z
		+ 4.166664568298827E-002) * z * z;
		*/
		p = coscof;
		y = *p++;
		y = y * z + *p++;
		y = y * z + *p++;
		y *= z * z;
		y -= 0.5 * z;
		y += 1.0;
	}
	else
	{
		/* Theoretical relative error = 3.8e-9 in [-pi/4, +pi/4] */
		/*
		y = ((-1.9515295891E-4 * z
		+ 8.3321608736E-3) * z
		- 1.6666654611E-1) * z * x;
		y += x;
		*/
		p = sincof;
		y = *p++;
		y = y * z + *p++;
		y = y * z + *p++;
		y *= z * x;
		y += x;
	}
	/*einitd();*/
	if ( sign < 0 )
		y = -y;
	return y;
}

extern "C" float cosf( float xx )
{
	float x, y, z;
	int j, sign;

	/* make argument positive */
	sign = 1;
	x = xx;
	if ( x < 0 )
		x = -x;

	if ( x > T24M1 )
	{
		// mtherr( "cosf", TLOSS );
		return( 0.0 );
	}

	j = FOPI * x; /* integer part of x/PIO4 */
	y = j;
	/* integer and fractional part modulo one octant */
	if ( j & 1 )	/* map zeros to origin */
	{
		j += 1;
		y += 1.0;
	}
	j &= 7;
	if ( j > 3 )
	{
		j -= 4;
		sign = -sign;
	}

	if ( j > 1 )
		sign = -sign;

	if ( x > lossth )
	{
		// mtherr( "cosf", PLOSS );
		x = x - y * PIO4F;
	}
	else
		/* Extended precision modular arithmetic */
		x = ( ( x - y * DP1 ) - y * DP2 ) - y * DP3;

	z = x * x;

	if ( ( j == 1 ) || ( j == 2 ) )
	{
		y = ( ( ( -1.9515295891E-4 * z
				  + 8.3321608736E-3 ) * z
				- 1.6666654611E-1 ) * z * x )
			+ x;
	}
	else
	{
		y = ( ( 2.443315711809948E-005 * z
				- 1.388731625493765E-003 ) * z
			  + 4.166664568298827E-002 ) * z * z;
		y -= 0.5 * z;
		y += 1.0;
	}
	if ( sign < 0 )
		y = -y;
	return y;
}

extern "C" float asinf( float xx )
{
	float a, x, z;
	int sign, flag;

	x = xx;

	if ( x > 0 )
	{
		sign = 1;
		a = x;
	}
	else
	{
		sign = -1;
		a = -x;
	}

	if ( a > 1.0 )
	{
		// mtherr( "asinf", DOMAIN );
		return( 0.0 );
	}

	if ( a < 1.0e-4 )
	{
		z = a;
		goto done;
	}

	if ( a > 0.5 )
	{
		z = 0.5 * ( 1.0 - a );
		x = sqrtf( z );
		flag = 1;
	}
	else
	{
		x = a;
		z = x * x;
		flag = 0;
	}

	z =
		( ( ( ( 4.2163199048E-2 * z
				+ 2.4181311049E-2 ) * z
			  + 4.5470025998E-2 ) * z
			+ 7.4953002686E-2 ) * z
		  + 1.6666752422E-1 ) * z * x
		+ x;

	if ( flag != 0 )
	{
		z = z + z;
		z = PIO2F - z;
	}
done:
	if ( sign < 0 )
		z = -z;
	return z;
}

extern "C" float acosf( float x )
{

	if ( x < -1.0 )
		goto domerr;

	if ( x < -0.5 )
		return( PIF - 2.0 * asinf( sqrtf( 0.5 * ( 1.0 + x ) ) ) );

	if ( x > 1.0 )
	{
		// domerr:	mtherr( "acosf", DOMAIN );
	domerr:
		return( 0.0 );
	}

	if ( x > 0.5 )
		return( 2.0 * asinf( sqrtf( 0.5 * ( 1.0 - x ) ) ) );

	return PIO2F - asinf( x );
}

extern "C" float atanf( float xx )
{
	float x, y, z;
	int sign;

	x = xx;

	/* make argument positive and save the sign */
	if ( xx < 0.0 )
	{
		sign = -1;
		x = -xx;
	}
	else
	{
		sign = 1;
		x = xx;
	}
	/* range reduction */
	if ( x > 2.414213562373095 )  /* tan 3pi/8 */
	{
		y = PIO2F;
		x = -( 1.0 / x );
	}

	else if ( x > 0.4142135623730950 ) /* tan pi/8 */
	{
		y = PIO4F;
		x = ( x - 1.0 ) / ( x + 1.0 );
	}
	else
		y = 0.0;

	z = x * x;
	y +=
		( ( ( 8.05374449538e-2 * z
			  - 1.38776856032E-1 ) * z
			+ 1.99777106478E-1 ) * z
		  - 3.33329491539E-1 ) * z * x
		+ x;

	if ( sign < 0 )
		y = -y;

	return( y );
}

extern "C" float atan2f( float y, float x )
{
	float z, w;
	int code;


	code = 0;

	if ( x < 0.0 )
		code = 2;
	if ( y < 0.0 )
		code |= 1;

	if ( x == 0.0 )
	{
		if ( code & 1 )
		{
			return( -PIO2F );
		}
		if ( y == 0.0 )
			return( 0.0 );
		return( PIO2F );
	}

	if ( y == 0.0 )
	{
		if ( code & 2 )
			return( PIF );
		return( 0.0 );
	}


	switch ( code )
	{
		default:
		case 0:
		case 1: w = 0.0; break;
		case 2: w = PIF; break;
		case 3: w = -PIF; break;
	}

	z = atanf( y / x );

	return( w + z );
}

bool isdigit( char ch )
{
	return ( ch >= '0' ) && ( ch <= '9' );
}

extern "C" double atof( const char* s )
{
	// This function stolen from either Rolf Neugebauer or Andrew Tolmach. 
	// Probably Rolf.
	double a = 0.0;
	int e = 0;
	int c;
	while ( ( c = *s++ ) != '\0' && isdigit( c ) )
	{
		a = a * 10.0 + ( c - '0' );
	}
	if ( c == '.' )
	{
		while ( ( c = *s++ ) != '\0' && isdigit( c ) )
		{
			a = a * 10.0 + ( c - '0' );
			e = e - 1;
		}
	}
	if ( c == 'e' || c == 'E' )
	{
		int sign = 1;
		int i = 0;
		c = *s++;
		if ( c == '+' )
			c = *s++;
		else if ( c == '-' )
		{
			c = *s++;
			sign = -1;
		}
		while ( isdigit( c ) )
		{
			i = i * 10 + ( c - '0' );
			c = *s++;
		}
		e += i * sign;
	}
	while ( e > 0 )
	{
		a *= 10.0;
		e--;
	}
	while ( e < 0 )
	{
		a *= 0.1;
		e++;
	}
	return a;
}

#define DBL_EPSILON 2.22044604925031308085e-16
static const double toint = 1 / DBL_EPSILON;
#define FORCE_EVAL(x) do {                        \
	if (sizeof(x) == sizeof(float)) {         \
		volatile float __x;               \
		__x = (x);                        \
	} else if (sizeof(x) == sizeof(double)) { \
		volatile double __x;              \
		__x = (x);                        \
	} else {                                  \
		volatile long double __x;         \
		__x = (x);                        \
	}                                         \
} while(0)

double ceil( double x )
{
	union
	{
		double f; uint64_t i;
	} u = { x };
	int e = u.i >> 52 & 0x7ff;
	double y;

	if ( e >= 0x3ff + 52 || x == 0 )
		return x;
	/* y = int(x) - x, where int(x) is an integer neighbor of x */
	if ( u.i >> 63 )
		y = x - toint + toint - x;
	else
		y = x + toint - toint - x;
	/* special case because of non-nearest rounding modes */
	if ( e <= 0x3ff - 1 )
	{
		FORCE_EVAL( y );
		return u.i >> 63 ? -0.0 : 1;
	}
	if ( y < 0 )
		return x + y + 1;
	return x + y;
}

extern "C" float ceilf( float x )
{
	return (float)ceil( x );
}

static __inline unsigned __FLOAT_BITS( float __f )
{
	union
	{
		float __f; unsigned __i;
	} __u;
	__u.__f = __f;
	return __u.__i;
}
static __inline unsigned long long __DOUBLE_BITS( double __f )
{
	union
	{
		double __f; unsigned long long __i;
	} __u;
	__u.__f = __f;
	return __u.__i;
}

#define FP_NAN       0
#define FP_INFINITE  1
#define FP_ZERO      2
#define FP_SUBNORMAL 3
#define FP_NORMAL    4

int __fpclassify( double x )
{
	union
	{
		double f; uint64_t i;
	} u = { x };
	int e = u.i >> 52 & 0x7ff;
	if ( !e ) return u.i << 1 ? FP_SUBNORMAL : FP_ZERO;
	if ( e == 0x7ff ) return u.i << 12 ? FP_NAN : FP_INFINITE;
	return FP_NORMAL;
}
int __fpclassifyl( long double x )
{
	return __fpclassify( x );
}

#define isnan(x) ( \
	sizeof(x) == sizeof(float) ? (__FLOAT_BITS(x) & 0x7fffffff) > 0x7f800000 : \
	sizeof(x) == sizeof(double) ? (__DOUBLE_BITS(x) & -1ULL>>1) > 0x7ffULL<<52 : \
	__fpclassifyl(x) == FP_NAN)


double fmod( double x, double y )
{
	union
	{
		double f; uint64_t i;
	} ux = { x }, uy = { y };
	int ex = ux.i >> 52 & 0x7ff;
	int ey = uy.i >> 52 & 0x7ff;
	int sx = ux.i >> 63;
	uint64_t i;

	/* in the followings uxi should be ux.i, but then gcc wrongly adds */
	/* float load/store to inner loops ruining performance and code size */
	uint64_t uxi = ux.i;

	if ( uy.i << 1 == 0 || isnan( y ) || ex == 0x7ff )
		return ( x * y ) / ( x * y );
	if ( uxi << 1 <= uy.i << 1 )
	{
		if ( uxi << 1 == uy.i << 1 )
			return 0 * x;
		return x;
	}

	/* normalize x and y */
	if ( !ex )
	{
		for ( i = uxi << 12; i >> 63 == 0; ex--, i <<= 1 );
		uxi <<= -ex + 1;
	}
	else
	{
		uxi &= -1ULL >> 12;
		uxi |= 1ULL << 52;
	}
	if ( !ey )
	{
		for ( i = uy.i << 12; i >> 63 == 0; ey--, i <<= 1 );
		uy.i <<= -ey + 1;
	}
	else
	{
		uy.i &= -1ULL >> 12;
		uy.i |= 1ULL << 52;
	}

	/* x mod y */
	for ( ; ex > ey; ex-- )
	{
		i = uxi - uy.i;
		if ( i >> 63 == 0 )
		{
			if ( i == 0 )
				return 0 * x;
			uxi = i;
		}
		uxi <<= 1;
	}
	i = uxi - uy.i;
	if ( i >> 63 == 0 )
	{
		if ( i == 0 )
			return 0 * x;
		uxi = i;
	}
	for ( ; uxi >> 52 == 0; uxi <<= 1, ex-- );

	/* scale result */
	if ( ex > 0 )
	{
		uxi -= 1ULL << 52;
		uxi |= (uint64_t)ex << 52;
	}
	else
	{
		uxi >>= -ex + 1;
	}
	uxi |= (uint64_t)sx << 63;
	ux.i = uxi;
	return ux.f;
}

extern "C" float fmodf( float x, float y )
{
	union
	{
		float f; uint32_t i;
	} ux = { x }, uy = { y };
	int ex = ux.i >> 23 & 0xff;
	int ey = uy.i >> 23 & 0xff;
	uint32_t sx = ux.i & 0x80000000;
	uint32_t i;
	uint32_t uxi = ux.i;

	if ( uy.i << 1 == 0 || isnan( y ) || ex == 0xff )
		return ( x * y ) / ( x * y );
	if ( uxi << 1 <= uy.i << 1 )
	{
		if ( uxi << 1 == uy.i << 1 )
			return 0 * x;
		return x;
	}

	/* normalize x and y */
	if ( !ex )
	{
		for ( i = uxi << 9; i >> 31 == 0; ex--, i <<= 1 );
		uxi <<= -ex + 1;
	}
	else
	{
		uxi &= -1U >> 9;
		uxi |= 1U << 23;
	}
	if ( !ey )
	{
		for ( i = uy.i << 9; i >> 31 == 0; ey--, i <<= 1 );
		uy.i <<= -ey + 1;
	}
	else
	{
		uy.i &= -1U >> 9;
		uy.i |= 1U << 23;
	}

	/* x mod y */
	for ( ; ex > ey; ex-- )
	{
		i = uxi - uy.i;
		if ( i >> 31 == 0 )
		{
			if ( i == 0 )
				return 0 * x;
			uxi = i;
		}
		uxi <<= 1;
	}
	i = uxi - uy.i;
	if ( i >> 31 == 0 )
	{
		if ( i == 0 )
			return 0 * x;
		uxi = i;
	}
	for ( ; uxi >> 23 == 0; uxi <<= 1, ex-- );

	/* scale result up */
	if ( ex > 0 )
	{
		uxi -= 1U << 23;
		uxi |= (uint32_t)ex << 23;
	}
	else
	{
		uxi >>= -ex + 1;
	}
	uxi |= sx;
	ux.i = uxi;
	return ux.f;
}

static const double
ln2_hi = 6.93147180369123816490e-01,  /* 3fe62e42 fee00000 */
ln2_lo = 1.90821492927058770002e-10,  /* 3dea39ef 35793c76 */
Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

extern "C" double log( double x )
{
	union
	{
		double f; uint64_t i;
	} u = { x };
	double hfsq, f, s, z, R, w, t1, t2, dk;
	uint32_t hx;
	int k;

	hx = u.i >> 32;
	k = 0;
	if ( hx < 0x00100000 || hx >> 31 )
	{
		if ( u.i << 1 == 0 )
			return -1 / ( x * x );  /* log(+-0)=-inf */
		if ( hx >> 31 )
			return ( x - x ) / 0.0; /* log(-#) = NaN */
							  /* subnormal number, scale x up */
		k -= 54;
		x *= 0x1p54;
		u.f = x;
		hx = u.i >> 32;
	}
	else if ( hx >= 0x7ff00000 )
	{
		return x;
	}
	else if ( hx == 0x3ff00000 && u.i << 32 == 0 )
		return 0;

	/* reduce x into [sqrt(2)/2, sqrt(2)] */
	hx += 0x3ff00000 - 0x3fe6a09e;
	k += (int)( hx >> 20 ) - 0x3ff;
	hx = ( hx & 0x000fffff ) + 0x3fe6a09e;
	u.i = (uint64_t)hx << 32 | ( u.i & 0xffffffff );
	x = u.f;

	f = x - 1.0;
	hfsq = 0.5 * f * f;
	s = f / ( 2.0 + f );
	z = s * s;
	w = z * z;
	t1 = w * ( Lg2 + w * ( Lg4 + w * Lg6 ) );
	t2 = z * ( Lg1 + w * ( Lg3 + w * ( Lg5 + w * Lg7 ) ) );
	R = t2 + t1;
	dk = k;
	return s * ( hfsq + R ) + dk * ln2_lo - hfsq + f + dk * ln2_hi;
}

static const double
bp[] = { 1.0, 1.5, },
dp_h[] = { 0.0, 5.84962487220764160156e-01, }, /* 0x3FE2B803, 0x40000000 */
dp_l[] = { 0.0, 1.35003920212974897128e-08, }, /* 0x3E4CFDEB, 0x43CFD006 */
two53 = 9007199254740992.0, /* 0x43400000, 0x00000000 */
huge = 1.0e300,
tiny = 1.0e-300,
/* poly coefs for (3/2)*(log(x)-2s-2/3*s**3 */
L1 = 5.99999999999994648725e-01, /* 0x3FE33333, 0x33333303 */
L2 = 4.28571428578550184252e-01, /* 0x3FDB6DB6, 0xDB6FABFF */
L3 = 3.33333329818377432918e-01, /* 0x3FD55555, 0x518F264D */
L4 = 2.72728123808534006489e-01, /* 0x3FD17460, 0xA91D4101 */
L5 = 2.30660745775561754067e-01, /* 0x3FCD864A, 0x93C9DB65 */
L6 = 2.06975017800338417784e-01, /* 0x3FCA7E28, 0x4A454EEF */
P1 = 1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
P2 = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
P3 = 6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
P4 = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
P5 = 4.13813679705723846039e-08, /* 0x3E663769, 0x72BEA4D0 */
lg2 = 6.93147180559945286227e-01, /* 0x3FE62E42, 0xFEFA39EF */
lg2_h = 6.93147182464599609375e-01, /* 0x3FE62E43, 0x00000000 */
lg2_l = -1.90465429995776804525e-09, /* 0xBE205C61, 0x0CA86C39 */
ovt = 8.0085662595372944372e-017, /* -(1024-log2(ovfl+.5ulp)) */
cp = 9.61796693925975554329e-01, /* 0x3FEEC709, 0xDC3A03FD =2/(3ln2) */
cp_h = 9.61796700954437255859e-01, /* 0x3FEEC709, 0xE0000000 =(float)cp */
cp_l = -7.02846165095275826516e-09, /* 0xBE3E2FE0, 0x145B01F5 =tail of cp_h*/
ivln2 = 1.44269504088896338700e+00, /* 0x3FF71547, 0x652B82FE =1/ln2 */
ivln2_h = 1.44269502162933349609e+00, /* 0x3FF71547, 0x60000000 =24b 1/ln2*/
ivln2_l = 1.92596299112661746887e-08; /* 0x3E54AE0B, 0xF85DDF44 =1/ln2 tail*/

#define EXTRACT_WORDS(hi,lo,d)                    \
do {                                              \
  union {double f; uint64_t i;} __u;              \
  __u.f = (d);                                    \
  (hi) = __u.i >> 32;                             \
  (lo) = (uint32_t)__u.i;                         \
} while (0)

#define SET_LOW_WORD(d,lo)                        \
do {                                              \
  union {double f; uint64_t i;} __u;              \
  __u.f = (d);                                    \
  __u.i &= 0xffffffff00000000ull;                 \
  __u.i |= (uint32_t)(lo);                        \
  (d) = __u.f;                                    \
} while (0)

#define SET_HIGH_WORD(d,hi)                       \
do {                                              \
  union {double f; uint64_t i;} __u;              \
  __u.f = (d);                                    \
  __u.i &= 0xffffffff;                            \
  __u.i |= (uint64_t)(hi) << 32;                  \
  (d) = __u.f;                                    \
} while (0)

#define GET_HIGH_WORD(hi,d)                       \
do {                                              \
  union {double f; uint64_t i;} __u;              \
  __u.f = (d);                                    \
  (hi) = __u.i >> 32;                             \
} while (0)

double scalbn( double x, int n )
{
	union
	{
		double f; uint64_t i;
	} u;
	double y = x;

	if ( n > 1023 )
	{
		y *= 0x1p1023;
		n -= 1023;
		if ( n > 1023 )
		{
			y *= 0x1p1023;
			n -= 1023;
			if ( n > 1023 )
				n = 1023;
		}
	}
	else if ( n < -1022 )
	{
		y *= 0x1p-1022;
		n += 1022;
		if ( n < -1022 )
		{
			y *= 0x1p-1022;
			n += 1022;
			if ( n < -1022 )
				n = -1022;
		}
	}
	u.i = (uint64_t)( 0x3ff + n ) << 52;
	x = y * u.f;
	return x;
}

#define INSERT_WORDS(d,hi,lo)                     \
do {                                              \
  union {double f; uint64_t i;} __u;              \
  __u.i = ((uint64_t)(hi)<<32) | (uint32_t)(lo);  \
  (d) = __u.f;                                    \
} while (0)


double sqrt( double x )
{
	double z;
	int32_t sign = (int)0x80000000;
	int32_t ix0, s0, q, m, t, i;
	uint32_t r, t1, s1, ix1, q1;

	EXTRACT_WORDS( ix0, ix1, x );

	/* take care of Inf and NaN */
	if ( ( ix0 & 0x7ff00000 ) == 0x7ff00000 )
	{
		return x * x + x;  /* sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN */
	}
	/* take care of zero */
	if ( ix0 <= 0 )
	{
		if ( ( ( ix0 & ~sign ) | ix1 ) == 0 )
			return x;  /* sqrt(+-0) = +-0 */
		if ( ix0 < 0 )
			return ( x - x ) / ( x - x );  /* sqrt(-ve) = sNaN */
	}
	/* normalize x */
	m = ix0 >> 20;
	if ( m == 0 )
	{  /* subnormal x */
		while ( ix0 == 0 )
		{
			m -= 21;
			ix0 |= ( ix1 >> 11 );
			ix1 <<= 21;
		}
		for ( i = 0; ( ix0 & 0x00100000 ) == 0; i++ )
			ix0 <<= 1;
		m -= i - 1;
		ix0 |= ix1 >> ( 32 - i );
		ix1 <<= i;
	}
	m -= 1023;    /* unbias exponent */
	ix0 = ( ix0 & 0x000fffff ) | 0x00100000;
	if ( m & 1 )
	{  /* odd m, double x to make it even */
		ix0 += ix0 + ( ( ix1 & sign ) >> 31 );
		ix1 += ix1;
	}
	m >>= 1;      /* m = [m/2] */

				  /* generate sqrt(x) bit by bit */
	ix0 += ix0 + ( ( ix1 & sign ) >> 31 );
	ix1 += ix1;
	q = q1 = s0 = s1 = 0;  /* [q,q1] = sqrt(x) */
	r = 0x00200000;        /* r = moving bit from right to left */

	while ( r != 0 )
	{
		t = s0 + r;
		if ( t <= ix0 )
		{
			s0 = t + r;
			ix0 -= t;
			q += r;
		}
		ix0 += ix0 + ( ( ix1 & sign ) >> 31 );
		ix1 += ix1;
		r >>= 1;
	}

	r = sign;
	while ( r != 0 )
	{
		t1 = s1 + r;
		t = s0;
		if ( t < ix0 || ( t == ix0 && t1 <= ix1 ) )
		{
			s1 = t1 + r;
			if ( ( t1 & sign ) == sign && ( s1 & sign ) == 0 )
				s0++;
			ix0 -= t;
			if ( ix1 < t1 )
				ix0--;
			ix1 -= t1;
			q1 += r;
		}
		ix0 += ix0 + ( ( ix1 & sign ) >> 31 );
		ix1 += ix1;
		r >>= 1;
	}

	/* use floating add to find out rounding direction */
	if ( ( ix0 | ix1 ) != 0 )
	{
		z = 1.0 - tiny; /* raise inexact flag */
		if ( z >= 1.0 )
		{
			z = 1.0 + tiny;
			if ( q1 == (uint32_t)0xffffffff )
			{
				q1 = 0;
				q++;
			}
			else if ( z > 1.0 )
			{
				if ( q1 == (uint32_t)0xfffffffe )
					q++;
				q1 += 2;
			}
			else
				q1 += q1 & 1;
		}
	}
	ix0 = ( q >> 1 ) + 0x3fe00000;
	ix1 = q1 >> 1;
	if ( q & 1 )
		ix1 |= sign;
	ix0 += m << 20;
	INSERT_WORDS( z, ix0, ix1 );
	return z;
}

double fabs( double x )
{
	union
	{
		double f; uint64_t i;
	} u = { x };
	u.i &= -1ULL / 2;
	return u.f;
}

extern "C" double pow( double x, double y )
{
	double z, ax, z_h, z_l, p_h, p_l;
	double y1, t1, t2, r, s, t, u, v, w;
	int32_t i, j, k, yisint, n;
	int32_t hx, hy, ix, iy;
	uint32_t lx, ly;

	EXTRACT_WORDS( hx, lx, x );
	EXTRACT_WORDS( hy, ly, y );
	ix = hx & 0x7fffffff;
	iy = hy & 0x7fffffff;

	/* x**0 = 1, even if x is NaN */
	if ( ( iy | ly ) == 0 )
		return 1.0;
	/* 1**y = 1, even if y is NaN */
	if ( hx == 0x3ff00000 && lx == 0 )
		return 1.0;
	/* NaN if either arg is NaN */
	if ( ix > 0x7ff00000 || ( ix == 0x7ff00000 && lx != 0 ) ||
		 iy > 0x7ff00000 || ( iy == 0x7ff00000 && ly != 0 ) )
		return x + y;

	/* determine if y is an odd int when x < 0
	* yisint = 0       ... y is not an integer
	* yisint = 1       ... y is an odd int
	* yisint = 2       ... y is an even int
	*/
	yisint = 0;
	if ( hx < 0 )
	{
		if ( iy >= 0x43400000 )
			yisint = 2; /* even integer y */
		else if ( iy >= 0x3ff00000 )
		{
			k = ( iy >> 20 ) - 0x3ff;  /* exponent */
			if ( k > 20 )
			{
				j = ly >> ( 52 - k );
				if ( ( j << ( 52 - k ) ) == ly )
					yisint = 2 - ( j & 1 );
			}
			else if ( ly == 0 )
			{
				j = iy >> ( 20 - k );
				if ( ( j << ( 20 - k ) ) == iy )
					yisint = 2 - ( j & 1 );
			}
		}
	}

	/* special value of y */
	if ( ly == 0 )
	{
		if ( iy == 0x7ff00000 )
		{  /* y is +-inf */
			if ( ( ( ix - 0x3ff00000 ) | lx ) == 0 )  /* (-1)**+-inf is 1 */
				return 1.0;
			else if ( ix >= 0x3ff00000 ) /* (|x|>1)**+-inf = inf,0 */
				return hy >= 0 ? y : 0.0;
			else                       /* (|x|<1)**+-inf = 0,inf */
				return hy >= 0 ? 0.0 : -y;
		}
		if ( iy == 0x3ff00000 )
		{    /* y is +-1 */
			if ( hy >= 0 )
				return x;
			y = 1 / x;
#if FLT_EVAL_METHOD!=0
			{
				union
				{
					double f; uint64_t i;
				} u = { y };
				uint64_t i = u.i & -1ULL / 2;
				if ( i >> 52 == 0 && ( i & ( i - 1 ) ) )
					FORCE_EVAL( (float)y );
			}
#endif
			return y;
		}
		if ( hy == 0x40000000 )    /* y is 2 */
			return x * x;
		if ( hy == 0x3fe00000 )
		{  /* y is 0.5 */
			if ( hx >= 0 )     /* x >= +0 */
				return sqrt( x );
		}
	}

	ax = fabs( x );
	/* special value of x */
	if ( lx == 0 )
	{
		if ( ix == 0x7ff00000 || ix == 0 || ix == 0x3ff00000 )
		{ /* x is +-0,+-inf,+-1 */
			z = ax;
			if ( hy < 0 )   /* z = (1/|x|) */
				z = 1.0 / z;
			if ( hx < 0 )
			{
				if ( ( ( ix - 0x3ff00000 ) | yisint ) == 0 )
				{
					z = ( z - z ) / ( z - z ); /* (-1)**non-int is NaN */
				}
				else if ( yisint == 1 )
					z = -z;          /* (x<0)**odd = -(|x|**odd) */
			}
			return z;
		}
	}

	s = 1.0; /* sign of result */
	if ( hx < 0 )
	{
		if ( yisint == 0 ) /* (x<0)**(non-int) is NaN */
			return ( x - x ) / ( x - x );
		if ( yisint == 1 ) /* (x<0)**(odd int) */
			s = -1.0;
	}

	/* |y| is huge */
	if ( iy > 0x41e00000 )
	{ /* if |y| > 2**31 */
		if ( iy > 0x43f00000 )
		{  /* if |y| > 2**64, must o/uflow */
			if ( ix <= 0x3fefffff )
				return hy < 0 ? huge * huge : tiny * tiny;
			if ( ix >= 0x3ff00000 )
				return hy > 0 ? huge * huge : tiny * tiny;
		}
		/* over/underflow if x is not close to one */
		if ( ix < 0x3fefffff )
			return hy < 0 ? s * huge * huge : s * tiny * tiny;
		if ( ix > 0x3ff00000 )
			return hy > 0 ? s * huge * huge : s * tiny * tiny;
		/* now |1-x| is tiny <= 2**-20, suffice to compute
		log(x) by x-x^2/2+x^3/3-x^4/4 */
		t = ax - 1.0;       /* t has 20 trailing zeros */
		w = ( t * t ) * ( 0.5 - t * ( 0.3333333333333333333333 - t * 0.25 ) );
		u = ivln2_h * t;      /* ivln2_h has 21 sig. bits */
		v = t * ivln2_l - w * ivln2;
		t1 = u + v;
		SET_LOW_WORD( t1, 0 );
		t2 = v - ( t1 - u );
	}
	else
	{
		double ss, s2, s_h, s_l, t_h, t_l;
		n = 0;
		/* take care subnormal number */
		if ( ix < 0x00100000 )
		{
			ax *= two53;
			n -= 53;
			GET_HIGH_WORD( ix, ax );
		}
		n += ( ( ix ) >> 20 ) - 0x3ff;
		j = ix & 0x000fffff;
		/* determine interval */
		ix = j | 0x3ff00000;   /* normalize ix */
		if ( j <= 0x3988E )      /* |x|<sqrt(3/2) */
			k = 0;
		else if ( j < 0xBB67A )  /* |x|<sqrt(3)   */
			k = 1;
		else
		{
			k = 0;
			n += 1;
			ix -= 0x00100000;
		}
		SET_HIGH_WORD( ax, ix );

		/* compute ss = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
		u = ax - bp[ k ];        /* bp[0]=1.0, bp[1]=1.5 */
		v = 1.0 / ( ax + bp[ k ] );
		ss = u * v;
		s_h = ss;
		SET_LOW_WORD( s_h, 0 );
		/* t_h=ax+bp[k] High */
		t_h = 0.0;
		SET_HIGH_WORD( t_h, ( ( ix >> 1 ) | 0x20000000 ) + 0x00080000 + ( k << 18 ) );
		t_l = ax - ( t_h - bp[ k ] );
		s_l = v * ( ( u - s_h * t_h ) - s_h * t_l );
		/* compute log(ax) */
		s2 = ss * ss;
		r = s2 * s2 * ( L1 + s2 * ( L2 + s2 * ( L3 + s2 * ( L4 + s2 * ( L5 + s2 * L6 ) ) ) ) );
		r += s_l * ( s_h + ss );
		s2 = s_h * s_h;
		t_h = 3.0 + s2 + r;
		SET_LOW_WORD( t_h, 0 );
		t_l = r - ( ( t_h - 3.0 ) - s2 );
		/* u+v = ss*(1+...) */
		u = s_h * t_h;
		v = s_l * t_h + t_l * ss;
		/* 2/(3log2)*(ss+...) */
		p_h = u + v;
		SET_LOW_WORD( p_h, 0 );
		p_l = v - ( p_h - u );
		z_h = cp_h * p_h;        /* cp_h+cp_l = 2/(3*log2) */
		z_l = cp_l * p_h + p_l * cp + dp_l[ k ];
		/* log2(ax) = (ss+..)*2/(3*log2) = n + dp_h + z_h + z_l */
		t = (double)n;
		t1 = ( ( z_h + z_l ) + dp_h[ k ] ) + t;
		SET_LOW_WORD( t1, 0 );
		t2 = z_l - ( ( ( t1 - t ) - dp_h[ k ] ) - z_h );
	}

	/* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
	y1 = y;
	SET_LOW_WORD( y1, 0 );
	p_l = ( y - y1 ) * t1 + y * t2;
	p_h = y1 * t1;
	z = p_l + p_h;
	EXTRACT_WORDS( j, i, z );
	if ( j >= 0x40900000 )
	{                      /* z >= 1024 */
		if ( ( ( j - 0x40900000 ) | i ) != 0 )        /* if z > 1024 */
			return s * huge * huge;         /* overflow */
		if ( p_l + ovt > z - p_h )
			return s * huge * huge;         /* overflow */
	}
	else if ( ( j & 0x7fffffff ) >= 0x4090cc00 )
	{  /* z <= -1075 */  // FIXME: instead of abs(j) use unsigned j
		if ( ( ( j - 0xc090cc00 ) | i ) != 0 )        /* z < -1075 */
			return s * tiny * tiny;         /* underflow */
		if ( p_l <= z - p_h )
			return s * tiny * tiny;         /* underflow */
	}
	/*
	* compute 2**(p_h+p_l)
	*/
	i = j & 0x7fffffff;
	k = ( i >> 20 ) - 0x3ff;
	n = 0;
	if ( i > 0x3fe00000 )
	{  /* if |z| > 0.5, set n = [z+0.5] */
		n = j + ( 0x00100000 >> ( k + 1 ) );
		k = ( ( n & 0x7fffffff ) >> 20 ) - 0x3ff;  /* new k for n */
		t = 0.0;
		SET_HIGH_WORD( t, n & ~( 0x000fffff >> k ) );
		n = ( ( n & 0x000fffff ) | 0x00100000 ) >> ( 20 - k );
		if ( j < 0 )
			n = -n;
		p_h -= t;
	}
	t = p_l + p_h;
	SET_LOW_WORD( t, 0 );
	u = t * lg2_h;
	v = ( p_l - ( t - p_h ) ) * lg2 + t * lg2_l;
	z = u + v;
	w = v - ( z - u );
	t = z * z;
	t1 = z - t * ( P1 + t * ( P2 + t * ( P3 + t * ( P4 + t * P5 ) ) ) );
	r = ( z * t1 ) / ( t1 - 2.0 ) - ( w + z * w );
	z = 1.0 - ( r - z );
	GET_HIGH_WORD( j, z );
	j += n << 20;
	if ( ( j >> 20 ) <= 0 )  /* subnormal output */
		z = scalbn( z, n );
	else
		SET_HIGH_WORD( z, j );
	return s * z;
}

#define GET_FLOAT_WORD(w,d)                       \
do {                                              \
  union {float f; uint32_t i;} __u;               \
  __u.f = (d);                                    \
  (w) = __u.i;                                    \
} while (0)

float fabsf( float x )
{
	union
	{
		float f; uint32_t i;
	} u = { x };
	u.i &= 0x7fffffff;
	return u.f;
}

#define SET_FLOAT_WORD(d,w)                       \
do {                                              \
  union {float f; uint32_t i;} __u;               \
  __u.i = (w);                                    \
  (d) = __u.f;                                    \
} while (0)

static const float two24 = 16777216.0;  /* 0x4b800000 */
float scalbnf( float x, int n )
{
	union
	{
		float f; uint32_t i;
	} u;
	float y = x;

	if ( n > 127 )
	{
		y *= 0x1p127f;
		n -= 127;
		if ( n > 127 )
		{
			y *= 0x1p127f;
			n -= 127;
			if ( n > 127 )
				n = 127;
		}
	}
	else if ( n < -126 )
	{
		y *= 0x1p-126f;
		n += 126;
		if ( n < -126 )
		{
			y *= 0x1p-126f;
			n += 126;
			if ( n < -126 )
				n = -126;
		}
	}
	u.i = (uint32_t)( 0x7f + n ) << 23;
	x = y * u.f;
	return x;
}

extern "C" float powf( float x, float y )
{
	float z, ax, z_h, z_l, p_h, p_l;
	float y1, t1, t2, r, s, sn, t, u, v, w;
	int32_t i, j, k, yisint, n;
	int32_t hx, hy, ix, iy, is;

	GET_FLOAT_WORD( hx, x );
	GET_FLOAT_WORD( hy, y );
	ix = hx & 0x7fffffff;
	iy = hy & 0x7fffffff;

	/* x**0 = 1, even if x is NaN */
	if ( iy == 0 )
		return 1.0f;
	/* 1**y = 1, even if y is NaN */
	if ( hx == 0x3f800000 )
		return 1.0f;
	/* NaN if either arg is NaN */
	if ( ix > 0x7f800000 || iy > 0x7f800000 )
		return x + y;

	/* determine if y is an odd int when x < 0
	* yisint = 0       ... y is not an integer
	* yisint = 1       ... y is an odd int
	* yisint = 2       ... y is an even int
	*/
	yisint = 0;
	if ( hx < 0 )
	{
		if ( iy >= 0x4b800000 )
			yisint = 2; /* even integer y */
		else if ( iy >= 0x3f800000 )
		{
			k = ( iy >> 23 ) - 0x7f;         /* exponent */
			j = iy >> ( 23 - k );
			if ( ( j << ( 23 - k ) ) == iy )
				yisint = 2 - ( j & 1 );
		}
	}

	/* special value of y */
	if ( iy == 0x7f800000 )
	{  /* y is +-inf */
		if ( ix == 0x3f800000 )      /* (-1)**+-inf is 1 */
			return 1.0f;
		else if ( ix > 0x3f800000 )  /* (|x|>1)**+-inf = inf,0 */
			return hy >= 0 ? y : 0.0f;
		else                       /* (|x|<1)**+-inf = 0,inf */
			return hy >= 0 ? 0.0f : -y;
	}
	if ( iy == 0x3f800000 )    /* y is +-1 */
		return hy >= 0 ? x : 1.0f / x;
	if ( hy == 0x40000000 )    /* y is 2 */
		return x * x;
	if ( hy == 0x3f000000 )
	{  /* y is  0.5 */
		if ( hx >= 0 )     /* x >= +0 */
			return sqrtf( x );
	}

	ax = fabsf( x );
	/* special value of x */
	if ( ix == 0x7f800000 || ix == 0 || ix == 0x3f800000 )
	{ /* x is +-0,+-inf,+-1 */
		z = ax;
		if ( hy < 0 )  /* z = (1/|x|) */
			z = 1.0f / z;
		if ( hx < 0 )
		{
			if ( ( ( ix - 0x3f800000 ) | yisint ) == 0 )
			{
				z = ( z - z ) / ( z - z ); /* (-1)**non-int is NaN */
			}
			else if ( yisint == 1 )
				z = -z;          /* (x<0)**odd = -(|x|**odd) */
		}
		return z;
	}

	sn = 1.0f; /* sign of result */
	if ( hx < 0 )
	{
		if ( yisint == 0 ) /* (x<0)**(non-int) is NaN */
			return ( x - x ) / ( x - x );
		if ( yisint == 1 ) /* (x<0)**(odd int) */
			sn = -1.0f;
	}

	/* |y| is huge */
	if ( iy > 0x4d000000 )
	{ /* if |y| > 2**27 */
	 /* over/underflow if x is not close to one */
		if ( ix < 0x3f7ffff8 )
			return hy < 0 ? sn * huge * huge : sn * tiny * tiny;
		if ( ix > 0x3f800007 )
			return hy > 0 ? sn * huge * huge : sn * tiny * tiny;
		/* now |1-x| is tiny <= 2**-20, suffice to compute
		log(x) by x-x^2/2+x^3/3-x^4/4 */
		t = ax - 1;     /* t has 20 trailing zeros */
		w = ( t * t ) * ( 0.5f - t * ( 0.333333333333f - t * 0.25f ) );
		u = ivln2_h * t;  /* ivln2_h has 16 sig. bits */
		v = t * ivln2_l - w * ivln2;
		t1 = u + v;
		GET_FLOAT_WORD( is, t1 );
		SET_FLOAT_WORD( t1, is & 0xfffff000 );
		t2 = v - ( t1 - u );
	}
	else
	{
		float s2, s_h, s_l, t_h, t_l;
		n = 0;
		/* take care subnormal number */
		if ( ix < 0x00800000 )
		{
			ax *= two24;
			n -= 24;
			GET_FLOAT_WORD( ix, ax );
		}
		n += ( ( ix ) >> 23 ) - 0x7f;
		j = ix & 0x007fffff;
		/* determine interval */
		ix = j | 0x3f800000;     /* normalize ix */
		if ( j <= 0x1cc471 )       /* |x|<sqrt(3/2) */
			k = 0;
		else if ( j < 0x5db3d7 )   /* |x|<sqrt(3)   */
			k = 1;
		else
		{
			k = 0;
			n += 1;
			ix -= 0x00800000;
		}
		SET_FLOAT_WORD( ax, ix );

		/* compute s = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
		u = ax - bp[ k ];   /* bp[0]=1.0, bp[1]=1.5 */
		v = 1.0f / ( ax + bp[ k ] );
		s = u * v;
		s_h = s;
		GET_FLOAT_WORD( is, s_h );
		SET_FLOAT_WORD( s_h, is & 0xfffff000 );
		/* t_h=ax+bp[k] High */
		is = ( ( ix >> 1 ) & 0xfffff000 ) | 0x20000000;
		SET_FLOAT_WORD( t_h, is + 0x00400000 + ( k << 21 ) );
		t_l = ax - ( t_h - bp[ k ] );
		s_l = v * ( ( u - s_h * t_h ) - s_h * t_l );
		/* compute log(ax) */
		s2 = s * s;
		r = s2 * s2 * ( L1 + s2 * ( L2 + s2 * ( L3 + s2 * ( L4 + s2 * ( L5 + s2 * L6 ) ) ) ) );
		r += s_l * ( s_h + s );
		s2 = s_h * s_h;
		t_h = 3.0f + s2 + r;
		GET_FLOAT_WORD( is, t_h );
		SET_FLOAT_WORD( t_h, is & 0xfffff000 );
		t_l = r - ( ( t_h - 3.0f ) - s2 );
		/* u+v = s*(1+...) */
		u = s_h * t_h;
		v = s_l * t_h + t_l * s;
		/* 2/(3log2)*(s+...) */
		p_h = u + v;
		GET_FLOAT_WORD( is, p_h );
		SET_FLOAT_WORD( p_h, is & 0xfffff000 );
		p_l = v - ( p_h - u );
		z_h = cp_h * p_h;  /* cp_h+cp_l = 2/(3*log2) */
		z_l = cp_l * p_h + p_l * cp + dp_l[ k ];
		/* log2(ax) = (s+..)*2/(3*log2) = n + dp_h + z_h + z_l */
		t = (float)n;
		t1 = ( ( ( z_h + z_l ) + dp_h[ k ] ) + t );
		GET_FLOAT_WORD( is, t1 );
		SET_FLOAT_WORD( t1, is & 0xfffff000 );
		t2 = z_l - ( ( ( t1 - t ) - dp_h[ k ] ) - z_h );
	}

	/* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
	GET_FLOAT_WORD( is, y );
	SET_FLOAT_WORD( y1, is & 0xfffff000 );
	p_l = ( y - y1 ) * t1 + y * t2;
	p_h = y1 * t1;
	z = p_l + p_h;
	GET_FLOAT_WORD( j, z );
	if ( j > 0x43000000 )          /* if z > 128 */
		return sn * huge * huge;  /* overflow */
	else if ( j == 0x43000000 )
	{  /* if z == 128 */
		if ( p_l + ovt > z - p_h )
			return sn * huge * huge;  /* overflow */
	}
	else if ( ( j & 0x7fffffff ) > 0x43160000 )  /* z < -150 */ // FIXME: check should be  (uint32_t)j > 0xc3160000
		return sn * tiny * tiny;  /* underflow */
	else if ( j == 0xc3160000 )
	{  /* z == -150 */
		if ( p_l <= z - p_h )
			return sn * tiny * tiny;  /* underflow */
	}
	/*
	* compute 2**(p_h+p_l)
	*/
	i = j & 0x7fffffff;
	k = ( i >> 23 ) - 0x7f;
	n = 0;
	if ( i > 0x3f000000 )
	{   /* if |z| > 0.5, set n = [z+0.5] */
		n = j + ( 0x00800000 >> ( k + 1 ) );
		k = ( ( n & 0x7fffffff ) >> 23 ) - 0x7f;  /* new k for n */
		SET_FLOAT_WORD( t, n & ~( 0x007fffff >> k ) );
		n = ( ( n & 0x007fffff ) | 0x00800000 ) >> ( 23 - k );
		if ( j < 0 )
			n = -n;
		p_h -= t;
	}
	t = p_l + p_h;
	GET_FLOAT_WORD( is, t );
	SET_FLOAT_WORD( t, is & 0xffff8000 );
	u = t * lg2_h;
	v = ( p_l - ( t - p_h ) ) * lg2 + t * lg2_l;
	z = u + v;
	w = v - ( z - u );
	t = z * z;
	t1 = z - t * ( P1 + t * ( P2 + t * ( P3 + t * ( P4 + t * P5 ) ) ) );
	r = ( z * t1 ) / ( t1 - 2.0f ) - ( w + z * w );
	z = 1.0f - ( r - z );
	GET_FLOAT_WORD( j, z );
	j += n << 23;
	if ( ( j >> 23 ) <= 0 )  /* subnormal output */
		z = scalbnf( z, n );
	else
		SET_FLOAT_WORD( z, j );
	return sn * z;
}

extern "C" float logf( float x )
{
	union
	{
		float f; uint32_t i;
	} u = { x };
	float hfsq, f, s, z, R, w, t1, t2, dk;
	uint32_t ix;
	int k;

	ix = u.i;
	k = 0;
	if ( ix < 0x00800000 || ix >> 31 )
	{  /* x < 2**-126  */
		if ( ix << 1 == 0 )
			return -1 / ( x * x );  /* log(+-0)=-inf */
		if ( ix >> 31 )
			return ( x - x ) / 0.0f; /* log(-#) = NaN */
							   /* subnormal number, scale up x */
		k -= 25;
		x *= 0x1p25f;
		u.f = x;
		ix = u.i;
	}
	else if ( ix >= 0x7f800000 )
	{
		return x;
	}
	else if ( ix == 0x3f800000 )
		return 0;

	/* reduce x into [sqrt(2)/2, sqrt(2)] */
	ix += 0x3f800000 - 0x3f3504f3;
	k += (int)( ix >> 23 ) - 0x7f;
	ix = ( ix & 0x007fffff ) + 0x3f3504f3;
	u.i = ix;
	x = u.f;

	f = x - 1.0f;
	s = f / ( 2.0f + f );
	z = s * s;
	w = z * z;
	t1 = w * ( Lg2 + w * Lg4 );
	t2 = z * ( Lg1 + w * Lg3 );
	R = t2 + t1;
	hfsq = 0.5f * f * f;
	dk = k;
	return s * ( hfsq + R ) + dk * ln2_lo - hfsq + f + dk * ln2_hi;
}

#define M_PI_2          1.57079632679489661923  /* pi/2 */
static const double T[] = {
	0x15554d3418c99f.0p-54, /* 0.333331395030791399758 */
	0x1112fd38999f72.0p-55, /* 0.133392002712976742718 */
	0x1b54c91d865afe.0p-57, /* 0.0533812378445670393523 */
	0x191df3908c33ce.0p-58, /* 0.0245283181166547278873 */
	0x185dadfcecf44e.0p-61, /* 0.00297435743359967304927 */
	0x1362b9bf971bcd.0p-59, /* 0.00946564784943673166728 */
};

float __tandf( double x, int odd )
{
	double z, r, w, s, t, u;

	z = x * x;
	/*
	* Split up the polynomial into small independent terms to give
	* opportunities for parallel evaluation.  The chosen splitting is
	* micro-optimized for Athlons (XP, X64).  It costs 2 multiplications
	* relative to Horner's method on sequential machines.
	*
	* We add the small terms from lowest degree up for efficiency on
	* non-sequential machines (the lowest degree terms tend to be ready
	* earlier).  Apart from this, we don't care about order of
	* operations, and don't need to to care since we have precision to
	* spare.  However, the chosen splitting is good for accuracy too,
	* and would give results as accurate as Horner's method if the
	* small terms were added from highest degree down.
	*/
	r = T[ 4 ] + z * T[ 5 ];
	t = T[ 2 ] + z * T[ 3 ];
	w = z * z;
	s = z * x;
	u = T[ 0 ] + z * T[ 1 ];
	r = ( x + s * u ) + ( s * w ) * ( t + w * r );
	return odd ? -1.0 / r : r;
}

extern "C" double floor(double x)
{
	union
	{
		double f; uint64_t i;
	} u = { x };
	int e = u.i >> 52 & 0x7ff;
	double y;

	if ( e >= 0x3ff + 52 || x == 0 )
		return x;
	/* y = int(x) - x, where int(x) is an integer neighbor of x */
	if ( u.i >> 63 )
		y = x - toint + toint - x;
	else
		y = x + toint - toint - x;
	/* special case because of non-nearest rounding modes */
	if ( e <= 0x3ff - 1 )
	{
		FORCE_EVAL( y );
		return u.i >> 63 ? -1 : 0;
	}
	if ( y > 0 )
		return x + y - 1;
	return x + y;
}

static const int init_jk[] = { 3,4,4,6 }; /* initial value for jk */
static const int32_t ipio2[] = {
	0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62,
	0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
	0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
	0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41,
	0x3991D6, 0x398353, 0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8,
	0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
	0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
	0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08,
	0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
	0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880,
	0x4D7327, 0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B,

#if LDBL_MAX_EXP > 1024
	0x47C419, 0xC367CD, 0xDCE809, 0x2A8359, 0xC4768B, 0x961CA6,
	0xDDAF44, 0xD15719, 0x053EA5, 0xFF0705, 0x3F7E33, 0xE832C2,
	0xDE4F98, 0x327DBB, 0xC33D26, 0xEF6B1E, 0x5EF89F, 0x3A1F35,
	0xCAF27F, 0x1D87F1, 0x21907C, 0x7C246A, 0xFA6ED5, 0x772D30,
	0x433B15, 0xC614B5, 0x9D19C3, 0xC2C4AD, 0x414D2C, 0x5D000C,
	0x467D86, 0x2D71E3, 0x9AC69B, 0x006233, 0x7CD2B4, 0x97A7B4,
	0xD55537, 0xF63ED7, 0x1810A3, 0xFC764D, 0x2A9D64, 0xABD770,
	0xF87C63, 0x57B07A, 0xE71517, 0x5649C0, 0xD9D63B, 0x3884A7,
	0xCB2324, 0x778AD6, 0x23545A, 0xB91F00, 0x1B0AF1, 0xDFCE19,
	0xFF319F, 0x6A1E66, 0x615799, 0x47FBAC, 0xD87F7E, 0xB76522,
	0x89E832, 0x60BFE6, 0xCDC4EF, 0x09366C, 0xD43F5D, 0xD7DE16,
	0xDE3B58, 0x929BDE, 0x2822D2, 0xE88628, 0x4D58E2, 0x32CAC6,
	0x16E308, 0xCB7DE0, 0x50C017, 0xA71DF3, 0x5BE018, 0x34132E,
	0x621283, 0x014883, 0x5B8EF5, 0x7FB0AD, 0xF2E91E, 0x434A48,
	0xD36710, 0xD8DDAA, 0x425FAE, 0xCE616A, 0xA4280A, 0xB499D3,
	0xF2A606, 0x7F775C, 0x83C2A3, 0x883C61, 0x78738A, 0x5A8CAF,
	0xBDD76F, 0x63A62D, 0xCBBFF4, 0xEF818D, 0x67C126, 0x45CA55,
	0x36D9CA, 0xD2A828, 0x8D61C2, 0x77C912, 0x142604, 0x9B4612,
	0xC459C4, 0x44C5C8, 0x91B24D, 0xF31700, 0xAD43D4, 0xE54929,
	0x10D5FD, 0xFCBE00, 0xCC941E, 0xEECE70, 0xF53E13, 0x80F1EC,
	0xC3E7B3, 0x28F8C7, 0x940593, 0x3E71C1, 0xB3092E, 0xF3450B,
	0x9C1288, 0x7B20AB, 0x9FB52E, 0xC29247, 0x2F327B, 0x6D550C,
	0x90A772, 0x1FE76B, 0x96CB31, 0x4A1679, 0xE27941, 0x89DFF4,
	0x9794E8, 0x84E6E2, 0x973199, 0x6BED88, 0x365F5F, 0x0EFDBB,
	0xB49A48, 0x6CA467, 0x427271, 0x325D8D, 0xB8159F, 0x09E5BC,
	0x25318D, 0x3974F7, 0x1C0530, 0x010C0D, 0x68084B, 0x58EE2C,
	0x90AA47, 0x02E774, 0x24D6BD, 0xA67DF7, 0x72486E, 0xEF169F,
	0xA6948E, 0xF691B4, 0x5153D1, 0xF20ACF, 0x339820, 0x7E4BF5,
	0x6863B2, 0x5F3EDD, 0x035D40, 0x7F8985, 0x295255, 0xC06437,
	0x10D86D, 0x324832, 0x754C5B, 0xD4714E, 0x6E5445, 0xC1090B,
	0x69F52A, 0xD56614, 0x9D0727, 0x50045D, 0xDB3BB4, 0xC576EA,
	0x17F987, 0x7D6B49, 0xBA271D, 0x296996, 0xACCCC6, 0x5414AD,
	0x6AE290, 0x89D988, 0x50722C, 0xBEA404, 0x940777, 0x7030F3,
	0x27FC00, 0xA871EA, 0x49C266, 0x3DE064, 0x83DD97, 0x973FA3,
	0xFD9443, 0x8C860D, 0xDE4131, 0x9D3992, 0x8C70DD, 0xE7B717,
	0x3BDF08, 0x2B3715, 0xA0805C, 0x93805A, 0x921110, 0xD8E80F,
	0xAF806C, 0x4BFFDB, 0x0F9038, 0x761859, 0x15A562, 0xBBCB61,
	0xB989C7, 0xBD4010, 0x04F2D2, 0x277549, 0xF6B6EB, 0xBB22DB,
	0xAA140A, 0x2F2689, 0x768364, 0x333B09, 0x1A940E, 0xAA3A51,
	0xC2A31D, 0xAEEDAF, 0x12265C, 0x4DC26D, 0x9C7A2D, 0x9756C0,
	0x833F03, 0xF6F009, 0x8C402B, 0x99316D, 0x07B439, 0x15200C,
	0x5BC3D8, 0xC492F5, 0x4BADC6, 0xA5CA4E, 0xCD37A7, 0x36A9E6,
	0x9492AB, 0x6842DD, 0xDE6319, 0xEF8C76, 0x528B68, 0x37DBFC,
	0xABA1AE, 0x3115DF, 0xA1AE00, 0xDAFB0C, 0x664D64, 0xB705ED,
	0x306529, 0xBF5657, 0x3AFF47, 0xB9F96A, 0xF3BE75, 0xDF9328,
	0x3080AB, 0xF68C66, 0x15CB04, 0x0622FA, 0x1DE4D9, 0xA4B33D,
	0x8F1B57, 0x09CD36, 0xE9424E, 0xA4BE13, 0xB52333, 0x1AAAF0,
	0xA8654F, 0xA5C1D2, 0x0F3F0B, 0xCD785B, 0x76F923, 0x048B7B,
	0x721789, 0x53A6C6, 0xE26E6F, 0x00EBEF, 0x584A9B, 0xB7DAC4,
	0xBA66AA, 0xCFCF76, 0x1D02D1, 0x2DF1B1, 0xC1998C, 0x77ADC3,
	0xDA4886, 0xA05DF7, 0xF480C6, 0x2FF0AC, 0x9AECDD, 0xBC5C3F,
	0x6DDED0, 0x1FC790, 0xB6DB2A, 0x3A25A3, 0x9AAF00, 0x9353AD,
	0x0457B6, 0xB42D29, 0x7E804B, 0xA707DA, 0x0EAA76, 0xA1597B,
	0x2A1216, 0x2DB7DC, 0xFDE5FA, 0xFEDB89, 0xFDBE89, 0x6C76E4,
	0xFCA906, 0x70803E, 0x156E85, 0xFF87FD, 0x073E28, 0x336761,
	0x86182A, 0xEABD4D, 0xAFE7B3, 0x6E6D8F, 0x396795, 0x5BBF31,
	0x48D784, 0x16DF30, 0x432DC7, 0x356125, 0xCE70C9, 0xB8CB30,
	0xFD6CBF, 0xA200A4, 0xE46C05, 0xA0DD5A, 0x476F21, 0xD21262,
	0x845CB9, 0x496170, 0xE0566B, 0x015299, 0x375550, 0xB7D51E,
	0xC4F133, 0x5F6E13, 0xE4305D, 0xA92E85, 0xC3B21D, 0x3632A1,
	0xA4B708, 0xD4B1EA, 0x21F716, 0xE4698F, 0x77FF27, 0x80030C,
	0x2D408D, 0xA0CD4F, 0x99A520, 0xD3A2B3, 0x0A5D2F, 0x42F9B4,
	0xCBDA11, 0xD0BE7D, 0xC1DB9B, 0xBD17AB, 0x81A2CA, 0x5C6A08,
	0x17552E, 0x550027, 0xF0147F, 0x8607E1, 0x640B14, 0x8D4196,
	0xDEBE87, 0x2AFDDA, 0xB6256B, 0x34897B, 0xFEF305, 0x9EBFB9,
	0x4F6A68, 0xA82A4A, 0x5AC44F, 0xBCF82D, 0x985AD7, 0x95C7F4,
	0x8D4D0D, 0xA63A20, 0x5F57A4, 0xB13F14, 0x953880, 0x0120CC,
	0x86DD71, 0xB6DEC9, 0xF560BF, 0x11654D, 0x6B0701, 0xACB08C,
	0xD0C0B2, 0x485551, 0x0EFB1E, 0xC37295, 0x3B06A3, 0x3540C0,
	0x7BDC06, 0xCC45E0, 0xFA294E, 0xC8CAD6, 0x41F3E8, 0xDE647C,
	0xD8649B, 0x31BED9, 0xC397A4, 0xD45877, 0xC5E369, 0x13DAF0,
	0x3C3ABA, 0x461846, 0x5F7555, 0xF5BDD2, 0xC6926E, 0x5D2EAC,
	0xED440E, 0x423E1C, 0x87C461, 0xE9FD29, 0xF3D6E7, 0xCA7C22,
	0x35916F, 0xC5E008, 0x8DD7FF, 0xE26A6E, 0xC6FDB0, 0xC10893,
	0x745D7C, 0xB2AD6B, 0x9D6ECD, 0x7B723E, 0x6A11C6, 0xA9CFF7,
	0xDF7329, 0xBAC9B5, 0x5100B7, 0x0DB2E2, 0x24BA74, 0x607DE5,
	0x8AD874, 0x2C150D, 0x0C1881, 0x94667E, 0x162901, 0x767A9F,
	0xBEFDFD, 0xEF4556, 0x367ED9, 0x13D9EC, 0xB9BA8B, 0xFC97C4,
	0x27A831, 0xC36EF1, 0x36C594, 0x56A8D8, 0xB5A8B4, 0x0ECCCF,
	0x2D8912, 0x34576F, 0x89562C, 0xE3CE99, 0xB920D6, 0xAA5E6B,
	0x9C2A3E, 0xCC5F11, 0x4A0BFD, 0xFBF4E1, 0x6D3B8E, 0x2C86E2,
	0x84D4E9, 0xA9B4FC, 0xD1EEEF, 0xC9352E, 0x61392F, 0x442138,
	0xC8D91B, 0x0AFC81, 0x6A4AFB, 0xD81C2F, 0x84B453, 0x8C994E,
	0xCC2254, 0xDC552A, 0xD6C6C0, 0x96190B, 0xB8701A, 0x649569,
	0x605A26, 0xEE523F, 0x0F117F, 0x11B5F4, 0xF5CBFC, 0x2DBC34,
	0xEEBC34, 0xCC5DE8, 0x605EDD, 0x9B8E67, 0xEF3392, 0xB817C9,
	0x9B5861, 0xBC57E1, 0xC68351, 0x103ED8, 0x4871DD, 0xDD1C2D,
	0xA118AF, 0x462C21, 0xD7F359, 0x987AD9, 0xC0549E, 0xFA864F,
	0xFC0656, 0xAE79E5, 0x362289, 0x22AD38, 0xDC9367, 0xAAE855,
	0x382682, 0x9BE7CA, 0xA40D51, 0xB13399, 0x0ED7A9, 0x480569,
	0xF0B265, 0xA7887F, 0x974C88, 0x36D1F9, 0xB39221, 0x4A827B,
	0x21CF98, 0xDC9F40, 0x5547DC, 0x3A74E1, 0x42EB67, 0xDF9DFE,
	0x5FD45E, 0xA4677B, 0x7AACBA, 0xA2F655, 0x23882B, 0x55BA41,
	0x086E59, 0x862A21, 0x834739, 0xE6E389, 0xD49EE5, 0x40FB49,
	0xE956FF, 0xCA0F1C, 0x8A59C5, 0x2BFA94, 0xC5C1D3, 0xCFC50F,
	0xAE5ADB, 0x86C547, 0x624385, 0x3B8621, 0x94792C, 0x876110,
	0x7B4C2A, 0x1A2C80, 0x12BF43, 0x902688, 0x893C78, 0xE4C4A8,
	0x7BDBE5, 0xC23AC4, 0xEAF426, 0x8A67F7, 0xBF920D, 0x2BA365,
	0xB1933D, 0x0B7CBD, 0xDC51A4, 0x63DD27, 0xDDE169, 0x19949A,
	0x9529A8, 0x28CE68, 0xB4ED09, 0x209F44, 0xCA984E, 0x638270,
	0x237C7E, 0x32B90F, 0x8EF5A7, 0xE75614, 0x08F121, 0x2A9DB5,
	0x4D7E6F, 0x5119A5, 0xABF9B5, 0xD6DF82, 0x61DD96, 0x023616,
	0x9F3AC4, 0xA1A283, 0x6DED72, 0x7A8D39, 0xA9B882, 0x5C326B,
	0x5B2746, 0xED3400, 0x7700D2, 0x55F4FC, 0x4D5901, 0x8071E0,
#endif
};

static const double PIo2[] = {
	1.57079625129699707031e+00, /* 0x3FF921FB, 0x40000000 */
	7.54978941586159635335e-08, /* 0x3E74442D, 0x00000000 */
	5.39030252995776476554e-15, /* 0x3CF84698, 0x80000000 */
	3.28200341580791294123e-22, /* 0x3B78CC51, 0x60000000 */
	1.27065575308067607349e-29, /* 0x39F01B83, 0x80000000 */
	1.22933308981111328932e-36, /* 0x387A2520, 0x40000000 */
	2.73370053816464559624e-44, /* 0x36E38222, 0x80000000 */
	2.16741683877804819444e-51, /* 0x3569F31D, 0x00000000 */
};

int __rem_pio2_large( double* x, double* y, int e0, int nx, int prec )
{
	int32_t jz, jx, jv, jp, jk, carry, n, iq[ 20 ], i, j, k, m, q0, ih;
	double z, fw, f[ 20 ], fq[ 20 ], q[ 20 ];

	/* initialize jk*/
	jk = init_jk[ prec ];
	jp = jk;

	/* determine jx,jv,q0, note that 3>q0 */
	jx = nx - 1;
	jv = ( e0 - 3 ) / 24;  if ( jv < 0 ) jv = 0;
	q0 = e0 - 24 * ( jv + 1 );

	/* set up f[0] to f[jx+jk] where f[jx+jk] = ipio2[jv+jk] */
	j = jv - jx; m = jx + jk;
	for ( i = 0; i <= m; i++, j++ )
		f[ i ] = j < 0 ? 0.0 : (double)ipio2[ j ];

	/* compute q[0],q[1],...q[jk] */
	for ( i = 0; i <= jk; i++ )
	{
		for ( j = 0, fw = 0.0; j <= jx; j++ )
			fw += x[ j ] * f[ jx + i - j ];
		q[ i ] = fw;
	}

	jz = jk;
recompute:
	/* distill q[] into iq[] reversingly */
	for ( i = 0, j = jz, z = q[ jz ]; j > 0; i++, j-- )
	{
		fw = (double)(int32_t)( 0x1p-24 * z );
		iq[ i ] = (int32_t)( z - 0x1p24 * fw );
		z = q[ j - 1 ] + fw;
	}

	/* compute n */
	z = scalbn( z, q0 );       /* actual value of z */
	z -= 8.0 * floor( z * 0.125 ); /* trim off integer >= 8 */
	n = (int32_t)z;
	z -= (double)n;
	ih = 0;
	if ( q0 > 0 )
	{  /* need iq[jz-1] to determine n */
		i = iq[ jz - 1 ] >> ( 24 - q0 ); n += i;
		iq[ jz - 1 ] -= i << ( 24 - q0 );
		ih = iq[ jz - 1 ] >> ( 23 - q0 );
	}
	else if ( q0 == 0 ) ih = iq[ jz - 1 ] >> 23;
	else if ( z >= 0.5 ) ih = 2;

	if ( ih > 0 )
	{  /* q > 0.5 */
		n += 1; carry = 0;
		for ( i = 0; i < jz; i++ )
		{  /* compute 1-q */
			j = iq[ i ];
			if ( carry == 0 )
			{
				if ( j != 0 )
				{
					carry = 1;
					iq[ i ] = 0x1000000 - j;
				}
			}
			else
				iq[ i ] = 0xffffff - j;
		}
		if ( q0 > 0 )
		{  /* rare case: chance is 1 in 12 */
			switch ( q0 )
			{
				case 1:
					iq[ jz - 1 ] &= 0x7fffff; break;
				case 2:
					iq[ jz - 1 ] &= 0x3fffff; break;
			}
		}
		if ( ih == 2 )
		{
			z = 1.0 - z;
			if ( carry != 0 )
				z -= scalbn( 1.0, q0 );
		}
	}

	/* check if recomputation is needed */
	if ( z == 0.0 )
	{
		j = 0;
		for ( i = jz - 1; i >= jk; i-- ) j |= iq[ i ];
		if ( j == 0 )
		{  /* need recomputation */
			for ( k = 1; iq[ jk - k ] == 0; k++ );  /* k = no. of terms needed */

			for ( i = jz + 1; i <= jz + k; i++ )
			{  /* add q[jz+1] to q[jz+k] */
				f[ jx + i ] = (double)ipio2[ jv + i ];
				for ( j = 0, fw = 0.0; j <= jx; j++ )
					fw += x[ j ] * f[ jx + i - j ];
				q[ i ] = fw;
			}
			jz += k;
			goto recompute;
		}
	}

	/* chop off zero terms */
	if ( z == 0.0 )
	{
		jz -= 1;
		q0 -= 24;
		while ( iq[ jz ] == 0 )
		{
			jz--;
			q0 -= 24;
		}
	}
	else
	{ /* break z into 24-bit if necessary */
		z = scalbn( z, -q0 );
		if ( z >= 0x1p24 )
		{
			fw = (double)(int32_t)( 0x1p-24 * z );
			iq[ jz ] = (int32_t)( z - 0x1p24 * fw );
			jz += 1;
			q0 += 24;
			iq[ jz ] = (int32_t)fw;
		}
		else
			iq[ jz ] = (int32_t)z;
	}

	/* convert integer "bit" chunk to floating-point value */
	fw = scalbn( 1.0, q0 );
	for ( i = jz; i >= 0; i-- )
	{
		q[ i ] = fw * (double)iq[ i ];
		fw *= 0x1p-24;
	}

	/* compute PIo2[0,...,jp]*q[jz,...,0] */
	for ( i = jz; i >= 0; i-- )
	{
		for ( fw = 0.0, k = 0; k <= jp && k <= jz - i; k++ )
			fw += PIo2[ k ] * q[ i + k ];
		fq[ jz - i ] = fw;
	}

	/* compress fq[] into y[] */
	switch ( prec )
	{
		case 0:
			fw = 0.0;
			for ( i = jz; i >= 0; i-- )
				fw += fq[ i ];
			y[ 0 ] = ih == 0 ? fw : -fw;
			break;
		case 1:
		case 2:
			fw = 0.0;
			for ( i = jz; i >= 0; i-- )
				fw += fq[ i ];
			// TODO: drop excess precision here once double_t is used
			fw = (double)fw;
			y[ 0 ] = ih == 0 ? fw : -fw;
			fw = fq[ 0 ] - fw;
			for ( i = 1; i <= jz; i++ )
				fw += fq[ i ];
			y[ 1 ] = ih == 0 ? fw : -fw;
			break;
		case 3:  /* painful */
			for ( i = jz; i > 0; i-- )
			{
				fw = fq[ i - 1 ] + fq[ i ];
				fq[ i ] += fq[ i - 1 ] - fw;
				fq[ i - 1 ] = fw;
			}
			for ( i = jz; i > 1; i-- )
			{
				fw = fq[ i - 1 ] + fq[ i ];
				fq[ i ] += fq[ i - 1 ] - fw;
				fq[ i - 1 ] = fw;
			}
			for ( fw = 0.0, i = jz; i >= 2; i-- )
				fw += fq[ i ];
			if ( ih == 0 )
			{
				y[ 0 ] = fq[ 0 ]; y[ 1 ] = fq[ 1 ]; y[ 2 ] = fw;
			}
			else
			{
				y[ 0 ] = -fq[ 0 ]; y[ 1 ] = -fq[ 1 ]; y[ 2 ] = -fw;
			}
	}
	return n & 7;
}

static const double
invpio2 = 6.36619772367581382433e-01, /* 0x3FE45F30, 0x6DC9C883 */
pio2_1 = 1.57079631090164184570e+00, /* 0x3FF921FB, 0x50000000 */
pio2_1t = 1.58932547735281966916e-08; /* 0x3E5110b4, 0x611A6263 */

int __rem_pio2f( float x, double* y )
{
	union
	{
		float f; uint32_t i;
	} u = { x };
	double tx[ 1 ], ty[ 1 ];
	double fn;
	uint32_t ix;
	int n, sign, e0;

	ix = u.i & 0x7fffffff;
	/* 25+53 bit pi is good enough for medium size */
	if ( ix < 0x4dc90fdb )
	{  /* |x| ~< 2^28*(pi/2), medium size */
	  /* Use a specialized rint() to get fn.  Assume round-to-nearest. */
		fn = x * invpio2 + toint - toint;
		n = (int32_t)fn;
		*y = x - fn * pio2_1 - fn * pio2_1t;
		return n;
	}
	if ( ix >= 0x7f800000 )
	{  /* x is inf or NaN */
		*y = x - x;
		return 0;
	}
	/* scale x into [2^23, 2^24-1] */
	sign = u.i >> 31;
	e0 = ( ix >> 23 ) - ( 0x7f + 23 );  /* e0 = ilogb(|x|)-23, positive */
	u.i = ix - ( e0 << 23 );
	tx[ 0 ] = u.f;
	n = __rem_pio2_large( tx, ty, e0, 1, 0 );
	if ( sign )
	{
		*y = -ty[ 0 ];
		return -n;
	}
	*y = ty[ 0 ];
	return n;
}

static const double
t1pio2 = 1 * M_PI_2, /* 0x3FF921FB, 0x54442D18 */
t2pio2 = 2 * M_PI_2, /* 0x400921FB, 0x54442D18 */
t3pio2 = 3 * M_PI_2, /* 0x4012D97C, 0x7F3321D2 */
t4pio2 = 4 * M_PI_2; /* 0x401921FB, 0x54442D18 */

extern "C" float tanf( float x )
{
	double y;
	uint32_t ix;
	unsigned n, sign;

	GET_FLOAT_WORD( ix, x );
	sign = ix >> 31;
	ix &= 0x7fffffff;

	if ( ix <= 0x3f490fda )
	{  /* |x| ~<= pi/4 */
		if ( ix < 0x39800000 )
		{  /* |x| < 2**-12 */
		  /* raise inexact if x!=0 and underflow if subnormal */
			FORCE_EVAL( ix < 0x00800000 ? x / 0x1p120f : x + 0x1p120f );
			return x;
		}
		return __tandf( x, 0 );
	}
	if ( ix <= 0x407b53d1 )
	{  /* |x| ~<= 5*pi/4 */
		if ( ix <= 0x4016cbe3 )  /* |x| ~<= 3pi/4 */
			return __tandf( ( sign ? x + t1pio2 : x - t1pio2 ), 1 );
		else
			return __tandf( ( sign ? x + t2pio2 : x - t2pio2 ), 0 );
	}
	if ( ix <= 0x40e231d5 )
	{  /* |x| ~<= 9*pi/4 */
		if ( ix <= 0x40afeddf )  /* |x| ~<= 7*pi/4 */
			return __tandf( ( sign ? x + t3pio2 : x - t3pio2 ), 1 );
		else
			return __tandf( ( sign ? x + t4pio2 : x - t4pio2 ), 0 );
	}

	/* tan(Inf or NaN) is NaN */
	if ( ix >= 0x7f800000 )
		return x - x;

	/* argument reduction */
	n = __rem_pio2f( x, &y );
	return __tandf( y, n & 1 );
}
