#include <windows.h>

extern "C" void abort()
{
	__debugbreak();
}

extern "C" EXCEPTION_DISPOSITION
__CxxFrameHandler3(
	IN PEXCEPTION_RECORD ExceptionRecord,
	IN PVOID EstablisherFrame,
	IN OUT PCONTEXT ContextRecord,
	IN OUT PDISPATCHER_CONTEXT DispatcherContext
)
{
	return {};
}

class exception;
class cxx_exception_type;
extern "C" void WINAPI _CxxThrowException( exception*, const cxx_exception_type* )
{}

extern "C" void __std_exception_copy( const struct __std_exception_data*, struct __std_exception_data* )
{}

extern "C" void __std_exception_destroy( struct __std_exception_data* )
{}

extern "C" void __std_terminate()
{}
