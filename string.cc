#include "printf/printf.hh"

extern "C" void abort();
extern "C" int atexit( [[maybe_unused]] void( *func )( ) );

extern "C" int __cdecl
strncmp(
	const char* lhs,
	const char* rhs,
	size_t count
)
{
	for ( ; count > 0; lhs++, rhs++, --count )
	{
		if ( *lhs != *rhs )
		{
			return ( ( *(unsigned char*)lhs < *(unsigned char*)rhs ) ? -1 : +1 );
		}
		else if ( *lhs == '\0' )
		{
			return 0;
		}
	}

	return 0;
}

extern "C" size_t strlen( const char* str )
{
	const char* char_ptr;
	const unsigned long int* longword_ptr;
	unsigned long int longword, himagic, lomagic;
	/* Handle the first few characters by reading one character at a time.
	Do this until CHAR_PTR is aligned on a longword boundary.  */
	for ( char_ptr = str; ( ( unsigned long int ) char_ptr
							& ( sizeof( longword ) - 1 ) ) != 0;
		  ++char_ptr )
		if ( *char_ptr == '\0' )
			return char_ptr - str;
	/* All these elucidatory comments refer to 4-byte longwords,
	but the theory applies equally well to 8-byte longwords.  */
	longword_ptr = ( unsigned long int* ) char_ptr;
	/* Bits 31, 24, 16, and 8 of this number are zero.  Call these bits
	the "holes."  Note that there is a hole just to the left of
	each byte, with an extra at the end:
	bits:  01111110 11111110 11111110 11111111
	bytes: AAAAAAAA BBBBBBBB CCCCCCCC DDDDDDDD
	The 1-bits make sure that carries propagate to the next 0-bit.
	The 0-bits provide holes for carries to fall into.  */
	himagic = 0x80808080L;
	lomagic = 0x01010101L;
	if ( sizeof( longword ) > 4 )
	{
		/* 64-bit version of the magic.  */
		/* Do the shift in two steps to avoid a warning if long has 32 bits.  */
		himagic = ( ( himagic << 16 ) << 16 ) | himagic;
		lomagic = ( ( lomagic << 16 ) << 16 ) | lomagic;
	}
	if ( sizeof( longword ) > 8 )
		abort();
	/* Instead of the traditional loop which tests each character,
	we will test a longword at a time.  The tricky part is testing
	if *any of the four* bytes in the longword in question are zero.  */
	for ( ;;)
	{
		longword = *longword_ptr++;
		if ( ( ( longword - lomagic ) & ~longword & himagic ) != 0 )
		{
			/* Which of the bytes was the zero?  If none of them were, it was
			a misfire; continue the search.  */
			const char* cp = (const char*)( longword_ptr - 1 );
			if ( cp[ 0 ] == 0 )
				return cp - str;
			if ( cp[ 1 ] == 0 )
				return cp - str + 1;
			if ( cp[ 2 ] == 0 )
				return cp - str + 2;
			if ( cp[ 3 ] == 0 )
				return cp - str + 3;
			if ( sizeof( longword ) > 4 )
			{
				if ( cp[ 4 ] == 0 )
					return cp - str + 4;
				if ( cp[ 5 ] == 0 )
					return cp - str + 5;
				if ( cp[ 6 ] == 0 )
					return cp - str + 6;
				if ( cp[ 7 ] == 0 )
					return cp - str + 7;
			}
		}
	}
}

extern "C" const char* __cdecl
strstr(
	const char* str,
	const char* substr
)
{
	char c;
	size_t len;

	c = *substr++;

	if ( !c )
	{
		return (char*)str;
	}

	len = strlen( substr );
	do
	{
		char sc;

		do
		{
			sc = *str++;

			if ( !sc )
			{
				return nullptr;
			}
		} while ( sc != c );
	} while ( strncmp( str, substr, len ) != 0 );

	return (char*)( str - 1 );
}

extern "C" int strcmp( const char* p1, const char* p2 )
{
	const unsigned char* s1 = (const unsigned char*)p1;
	const unsigned char* s2 = (const unsigned char*)p2;
	unsigned char c1, c2;
	do
	{
		c1 = (unsigned char)*s1++;
		c2 = (unsigned char)*s2++;
		if ( c1 == '\0' )
			return c1 - c2;
	} while ( c1 == c2 );
	return c1 - c2;
}

extern "C" size_t
wcslen( const wchar_t* s )
{
	size_t len = 0;
	while ( s[ len ] != L'\0' )
	{
		if ( s[ ++len ] == L'\0' )
			return len;
		if ( s[ ++len ] == L'\0' )
			return len;
		if ( s[ ++len ] == L'\0' )
			return len;
		++len;
	}
	return len;
}


extern "C" int isupper( int c )
{
	return ( c >= 'A' && c <= 'Z' );
}

extern "C" int tolower( int c )
{
	return isupper( c ) ? (c)-'A' + 'a' : c;
}

extern "C" unsigned int towlower( unsigned int c )
{
	return ( c < 0x00ff ? (unsigned int)( tolower( (int)c ) ) : c );
}

