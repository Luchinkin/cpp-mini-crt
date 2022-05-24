#include <Windows.h>
#include <assert.h>
#include <new>
#include <stddef.h>
#include <stdlib.h>

//#define malloc(x) HeapAlloc(GetProcessHeap(), 0, x)
//#define realloc(x,s) HeapReAlloc(GetProcessHeap(), 0, x, s)
//#define free(x)     HeapFree(GetProcessHeap(), 0, x)

#ifdef _USE_WIN_HEAP

void* malloc( size_t sizemem )
{
	return HeapAlloc( GetProcessHeap(), 0, sizemem );
}

void* realloc( void* ptr, size_t newsize )
{
	return HeapReAlloc( GetProcessHeap(), 0, ptr, newsize );
}

void free( void* ptr )
{
	HeapFree( GetProcessHeap(), 0, ptr );
}

#else


extern "C" void* malloc( size_t sizemem );
extern "C" void* realloc( void* ptr, size_t newsize );
extern "C" void free( void* ptr );

#endif

void* operator new( size_t size )
{
	void* const p = ::operator new( size, std::nothrow );

	return p;
}

void* operator new( size_t, std::align_val_t )
{
	abort();
}

void* operator new( size_t size, const std::nothrow_t& ) noexcept
{
	if ( size == 0 )
	{
		++size;
	}
	return malloc( size );
}

void* operator new( size_t, std::align_val_t, const std::nothrow_t& ) noexcept
{
	abort();
}

void operator delete( void* ptr ) noexcept
{
	::operator delete( ptr, std::nothrow );
}

void operator delete( void* ptr, size_t ) noexcept
{
	::operator delete( ptr, std::nothrow );
}

void operator delete( void*, std::align_val_t ) noexcept
{
	abort();
}

void operator delete( void* ptr, const std::nothrow_t& ) noexcept
{
	if ( ptr )
	{
		free( ptr );
	}
}

void operator delete( void*, std::align_val_t, const std::nothrow_t& ) noexcept
{
	abort();
}

void* operator new[]( size_t size )
{
	return ::operator new( size );
}

void* operator new[]( size_t, std::align_val_t )
{
	abort();
}

void* operator new[]( size_t size, const std::nothrow_t& ) noexcept
{
	return ::operator new( size, std::nothrow );
}

void* operator new[]( size_t, std::align_val_t, const std::nothrow_t& ) noexcept
{
	abort();
}

void operator delete[]( void* ptr ) noexcept
{
	::operator delete( ptr );
}

void operator delete[]( void* ptr, size_t size ) noexcept
{
	::operator delete( ptr, size );
}

void operator delete[]( void*, std::align_val_t ) noexcept
{
	abort();
}

void operator delete[]( void*, size_t, std::align_val_t ) noexcept
{
	abort();
}

void operator delete[]( void* ptr, const std::nothrow_t& ) noexcept
{
	::operator delete( ptr, std::nothrow );
}

void operator delete[]( void*, std::align_val_t, const std::nothrow_t& ) noexcept
{
	abort();
}

int __cdecl
memcmp(
	const void* lhs,
	const void* rhs,
	size_t count
)
{
	unsigned char u1, u2;
	unsigned char* s1 = (unsigned char*)lhs;
	unsigned char* s2 = (unsigned char*)rhs;

	for ( ; count--; s1++, s2++ )
	{
		u1 = *(unsigned char*)s1;
		u2 = *(unsigned char*)s2;

		if ( u1 != u2 )
		{
			return ( u1 - u2 );
		}
	}

	return 0;
}

void* __cdecl
memset(
	void* dest,
	int ch,
	size_t count
)
{
	if ( count )
	{
		char* d = (char*)dest;

		do
		{
			*d++ = ch;
		} while ( --count );
	}

	return dest;
}

void* __cdecl
memcpy(
	void* dest,
	const void* src,
	size_t count
)
{
	char* d = (char*)dest;
	char* s = (char*)src;

	while ( count-- )
	{
		*d++ = *s++;
	}

	return dest;
}

void* __cdecl
memmove(
	void* dest,
	const void* src,
	size_t count
)
{
	const char* s = (const char*)src;
	char* d = (char*)dest;

	if ( !count )
	{
		return dest;
	}

	if ( dest <= s )
	{
		return memcpy( dest, s, count );
	}

	s += count;
	d += count;

	while ( count-- )
	{
		*--d = *--s;
	}

	return dest;
}
