#include <test_check.hpp>

#include <iostream>

int main()
{
    MALLOY_CHECK_TRUE(true);
    MALLOY_CHECK_EQ(1 + 1, 2);

    std::cout << "malloy_smoke_tests passed\n";
    return 0;
}
