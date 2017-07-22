#pragma once

namespace chm
{

#define UnitTest(x) \
	do { \
		if (!(x)) { \
			std::cout << "TEST_FAILED: " << #x << std::endl; \
			std::abort(); \
		} \
	} while (0)

namespace unit_tests
{
	void testAll();
}

}
