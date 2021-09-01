#include <cstdint>

extern "C" void _invalid_parameter_noinfo_noreturn()
{
}

extern "C" void __security_init_cookie()
{
}

extern "C" void __security_check_cookie( uintptr_t StackCookie )
{
}

extern "C" uintptr_t __security_cookie = 0x0;

