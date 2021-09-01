#include <stddef.h>
#include <stdlib.h>
#include <exception>

// #undef _HAS_EXCEPTIONS
// #define _HAS_EXCEPTIONS 0

inline void RAISE_EXCEPTION()
{
	__debugbreak();
}

namespace std
{
	void __cdecl
		_Xlength_error( _In_z_ const char* )
	{
		RAISE_EXCEPTION();
	}

	void __cdecl
		_Xbad_alloc()
	{
		RAISE_EXCEPTION();
	}

	void __cdecl
		_Xinvalid_argument( _In_z_ const char* )
	{
		RAISE_EXCEPTION();
	}

	void __cdecl
		_Xout_of_range( _In_z_ const char* )
	{
		RAISE_EXCEPTION();
	}

	void __cdecl
		_Xoverflow_error( _In_z_ const char* )
	{
		RAISE_EXCEPTION();
	}

	void __cdecl
		_Xruntime_error( _In_z_ const char* )
	{
		RAISE_EXCEPTION();
	}


	void __cdecl
		my_raise_handler( class stdext::exception const& )
	{
		RAISE_EXCEPTION();
	}
	
	_Prhand _Raise_handler = my_raise_handler;

	void __cdecl
		_Xbad_function_call( void )
	{
		RAISE_EXCEPTION();
	}


}

extern "C" int __cdecl atexit(
	[[maybe_unused]] void( __cdecl * func )( void )
)
{
	return 0;
}

extern "C" void __fastcall __GSHandlerCheck()
{}
extern "C" void __fastcall __guard_dispatch_icall_fptr()
{}
extern "C" void __fastcall __C_specific_handler()
{}
extern "C" void __fastcall __isa_available()
{}
extern "C" void __fastcall _local_unwind()
{}
extern "C" void __fastcall __report_rangecheckfailure()
{}

void __cdecl
_invoke_watson(
	[[maybe_unused]] wchar_t const* const expression,
	[[maybe_unused]] wchar_t const* const function_name,
	[[maybe_unused]] wchar_t const* const file_name,
	[[maybe_unused]] unsigned int   const line_number,
	[[maybe_unused]] uintptr_t      const reserved )
{
	RAISE_EXCEPTION();
}
