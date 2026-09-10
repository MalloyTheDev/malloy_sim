#pragma once

#include <cmath>
#include <iostream>

#define MALLOY_CHECK_TRUE(expression)                                                  \
    do                                                                                  \
    {                                                                                   \
        if (!(expression))                                                              \
        {                                                                               \
            std::cerr << "MALLOY_CHECK_TRUE failed: " #expression                       \
                      << " at " << __FILE__ << ":" << __LINE__ << '\n';                 \
            return 1;                                                                   \
        }                                                                               \
    } while (false)

#define MALLOY_CHECK_FALSE(expression) MALLOY_CHECK_TRUE(!(expression))

#define MALLOY_CHECK_EQ(left, right)                                                    \
    do                                                                                  \
    {                                                                                   \
        const auto malloy_left_value = (left);                                          \
        const auto malloy_right_value = (right);                                        \
        if (!(malloy_left_value == malloy_right_value))                                 \
        {                                                                               \
            std::cerr << "MALLOY_CHECK_EQ failed: " #left " == " #right                 \
                      << " at " << __FILE__ << ":" << __LINE__                          \
                      << " left=" << malloy_left_value                                  \
                      << " right=" << malloy_right_value << '\n';                       \
            return 1;                                                                   \
        }                                                                               \
    } while (false)

// Fails when the two values differ by more than epsilon, and ALSO when either
// one is NaN. The comparison is negated deliberately: every comparison
// involving a NaN is false, so the direct form (difference > epsilon) is false
// for a NaN and would silently pass. !(difference <= epsilon) is true for it.
#define MALLOY_CHECK_NEAR(left, right, epsilon)                                         \
    do                                                                                  \
    {                                                                                   \
        const auto malloy_left_value = (left);                                          \
        const auto malloy_right_value = (right);                                        \
        const auto malloy_epsilon_value = (epsilon);                                    \
        if (!(std::abs(malloy_left_value - malloy_right_value) <=                       \
              malloy_epsilon_value))                                                    \
        {                                                                               \
            std::cerr << "MALLOY_CHECK_NEAR failed: " #left " ~= " #right               \
                      << " at " << __FILE__ << ":" << __LINE__                          \
                      << " left=" << malloy_left_value                                  \
                      << " right=" << malloy_right_value                                \
                      << " epsilon=" << malloy_epsilon_value << '\n';                   \
            return 1;                                                                   \
        }                                                                               \
    } while (false)

// Component-wise near-equality for any type exposing .x and .y members.
// Deliberately does not include vec2.hpp, so this generic test header stays
// decoupled from the math module; it just expands where such a type is in scope.
// Negated per component, for the NaN reason given above.
#define MALLOY_CHECK_VEC2_NEAR(left, right, epsilon)                                    \
    do                                                                                  \
    {                                                                                   \
        const auto malloy_vec2_left = (left);                                           \
        const auto malloy_vec2_right = (right);                                         \
        const auto malloy_vec2_epsilon = (epsilon);                                     \
        if (!(std::abs(malloy_vec2_left.x - malloy_vec2_right.x) <=                     \
              malloy_vec2_epsilon) ||                                                   \
            !(std::abs(malloy_vec2_left.y - malloy_vec2_right.y) <=                     \
              malloy_vec2_epsilon))                                                     \
        {                                                                               \
            std::cerr << "MALLOY_CHECK_VEC2_NEAR failed: " #left " ~= " #right          \
                      << " at " << __FILE__ << ":" << __LINE__                          \
                      << " left=(" << malloy_vec2_left.x << ", "                        \
                      << malloy_vec2_left.y << ")"                                      \
                      << " right=(" << malloy_vec2_right.x << ", "                      \
                      << malloy_vec2_right.y << ")"                                      \
                      << " epsilon=" << malloy_vec2_epsilon << '\n';                    \
            return 1;                                                                   \
        }                                                                               \
    } while (false)
