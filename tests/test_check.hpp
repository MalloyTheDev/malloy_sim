#pragma once

#include <cmath>
#include <iostream>

#define MALLOY_CHECK_TRUE(expression)                                                        \
    do                                                                                        \
    {                                                                                         \
        if (!(expression))                                                                    \
        {                                                                                     \
            std::cerr << "MALLOY_CHECK_TRUE failed: " #expression                            \
                      << " at " << __FILE__ << ":" << __LINE__ << '\n';                     \
            return 1;                                                                         \
        }                                                                                     \
    } while (false)

#define MALLOY_CHECK_FALSE(expression) MALLOY_CHECK_TRUE(!(expression))

#define MALLOY_CHECK_EQ(left, right)                                                          \
    do                                                                                        \
    {                                                                                         \
        const auto malloy_left_value = (left);                                                \
        const auto malloy_right_value = (right);                                              \
        if (!(malloy_left_value == malloy_right_value))                                       \
        {                                                                                     \
            std::cerr << "MALLOY_CHECK_EQ failed: " #left " == " #right                     \
                      << " at " << __FILE__ << ":" << __LINE__                              \
                      << " left=" << malloy_left_value                                       \
                      << " right=" << malloy_right_value << '\n';                            \
            return 1;                                                                         \
        }                                                                                     \
    } while (false)

#define MALLOY_CHECK_NEAR(left, right, epsilon)                                                \
    do                                                                                        \
    {                                                                                         \
        const auto malloy_left_value = (left);                                                \
        const auto malloy_right_value = (right);                                              \
        const auto malloy_epsilon_value = (epsilon);                                          \
        if (std::abs(malloy_left_value - malloy_right_value) > malloy_epsilon_value)          \
        {                                                                                     \
            std::cerr << "MALLOY_CHECK_NEAR failed: " #left " ~= " #right                   \
                      << " at " << __FILE__ << ":" << __LINE__                              \
                      << " left=" << malloy_left_value                                       \
                      << " right=" << malloy_right_value                                     \
                      << " epsilon=" << malloy_epsilon_value << '\n';                        \
            return 1;                                                                         \
        }                                                                                     \
    } while (false)
