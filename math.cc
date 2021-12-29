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

extern "C" double log(double x)
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

extern "C" float logf(float x)
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
