#pragma once

#include <iostream>

/* empty "real" statement */
#define CHC_NOP do { } while (0)

#ifdef NDEBUG
#define Debug(x) CHC_NOP
#else
#define Debug(x) do { std::cout << x << std::endl; } while (0)
#endif

#if defined(NVERBOSE) && defined(NDEBUG)
#define Print(x) CHC_NOP
#else
#define Print(x) do { std::cout << x << std::endl; } while (0)
#endif

#if defined(NVERBOSE) && defined(NDEBUG)
#define WorkingSign(x) CHC_NOP
#else
#define WorkingSign() do { std::cout << "+" << std::flush; } while (0)
#endif

#define Unused(x) ((void)x)

#ifdef NDEBUG
# undef NDEBUG
# include <cassert>
# define NDEBUG
# define debug_assert(...) CHC_NOP
#else
# include <cassert>
# define debug_assert(...) assert(__VA_ARGS__)
#endif

namespace chm
{
	typedef unsigned int uint;
}
