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
static float PIO2F =  1.570796326794896619F;

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

	if( x > 0 )
	{
		sign = 1;
		a = x;
	}
	else
	{
		sign = -1;
		a = -x;
	}

	if( a > 1.0 )
	{
		// mtherr( "asinf", DOMAIN );
		return( 0.0 );
	}

	if( a < 1.0e-4 )
	{
		z = a;
		goto done;
	}

	if( a > 0.5 )
	{
		z = 0.5 * (1.0 - a);
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
		(((( 4.2163199048E-2 * z
			 + 2.4181311049E-2) * z
			+ 4.5470025998E-2) * z
		   + 7.4953002686E-2) * z
		  + 1.6666752422E-1) * z * x
		+ x;

	if( flag != 0 )
	{
		z = z + z;
		z = PIO2F - z;
	}
done:
	if( sign < 0 )
		z = -z;
	return z;
}

extern "C" float acosf( float x )
{

	if( x < -1.0 )
		goto domerr;

	if( x < -0.5) 
		return( PIF - 2.0 * asinf( sqrtf(0.5*(1.0+x)) ) );

	if( x > 1.0 )
	{
	// domerr:	mtherr( "acosf", DOMAIN );
		domerr:
		return( 0.0 );
	}

	if( x > 0.5 )
		return( 2.0 * asinf(  sqrtf(0.5*(1.0-x) ) ) );

	return PIO2F - asinf(x);
}

extern "C" float atanf( float xx )
{
	float x, y, z;
	int sign;

	x = xx;

	/* make argument positive and save the sign */
	if( xx < 0.0 )
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
	if( x > 2.414213562373095 )  /* tan 3pi/8 */
	{
		y = PIO2F;
		x = -( 1.0/x );
	}

	else if( x > 0.4142135623730950 ) /* tan pi/8 */
	{
		y = PIO4F;
		x = (x-1.0)/(x+1.0);
	}
	else
		y = 0.0;

	z = x * x;
	y +=
		((( 8.05374449538e-2 * z
			- 1.38776856032E-1) * z
		   + 1.99777106478E-1) * z
		  - 3.33329491539E-1) * z * x
		+ x;

	if( sign < 0 )
		y = -y;

	return( y );
}

extern "C" float atan2f( float y, float x )
{
	float z, w;
	int code;


	code = 0;

	if( x < 0.0 )
		code = 2;
	if( y < 0.0 )
		code |= 1;

	if( x == 0.0 )
	{
		if( code & 1 )
		{
			return( -PIO2F );
		}
		if( y == 0.0 )
			return( 0.0 );
		return( PIO2F );
	}

	if( y == 0.0 )
	{
		if( code & 2 )
			return( PIF );
		return( 0.0 );
	}


	switch( code )
	{
		default:
		case 0:
		case 1: w = 0.0; break;
		case 2: w = PIF; break;
		case 3: w = -PIF; break;
	}

	z = atanf( y/x );

	return( w + z );
}

