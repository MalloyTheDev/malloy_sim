#include <test_check.hpp>

#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace
{

// Redirects std::cerr for as long as it is alive. A check that is SUPPOSED to
// fail still reports itself, and that report would otherwise look like a real
// failure in the test log.
class CapturedCerr
{
public:
    CapturedCerr() : previous_{std::cerr.rdbuf(sink_.rdbuf())} {}
    ~CapturedCerr() { std::cerr.rdbuf(previous_); }

    CapturedCerr(const CapturedCerr&) = delete;
    CapturedCerr& operator=(const CapturedCerr&) = delete;
    CapturedCerr(CapturedCerr&&) = delete;
    CapturedCerr& operator=(CapturedCerr&&) = delete;

    std::string text() const { return sink_.str(); }

private:
    std::ostringstream sink_;
    std::streambuf* previous_;
};

// A minimal type with .x and .y. MALLOY_CHECK_VEC2_NEAR deliberately does not
// know about malloy::math::Vec2, and this executable deliberately does not link
// it, so the macro is exercised exactly as decoupled as it claims to be.
struct XY
{
    double x{};
    double y{};
};

// The macros return 1 from the enclosing function when they fail, so a failing
// case has to run inside its own function and be judged by its return value.
int near_case(double left, double right, double epsilon)
{
    MALLOY_CHECK_NEAR(left, right, epsilon);
    return 0;
}

int vec2_near_case(XY left, XY right, double epsilon)
{
    MALLOY_CHECK_VEC2_NEAR(left, right, epsilon);
    return 0;
}

// std::cerr is captured around the case itself rather than around the
// assertions in main, so that a genuine failure of THOSE is still printed.
int quiet_near(double left, double right, double epsilon,
               std::string* message = nullptr)
{
    CapturedCerr captured;
    const int result = near_case(left, right, epsilon);
    if (message != nullptr)
    {
        *message = captured.text();
    }
    return result;
}

int quiet_vec2_near(XY left, XY right, double epsilon)
{
    CapturedCerr captured;
    return vec2_near_case(left, right, epsilon);
}

} // namespace

int main()
{
    MALLOY_CHECK_TRUE(true);
    MALLOY_CHECK_EQ(1 + 1, 2);

    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();

    // --- The harness itself. Every other test executable in this project
    //     trusts these macros, so a hole in one of them is invisible in all of
    //     them at once. ---
    {
        // Ordinary behaviour: inside the tolerance passes, outside fails, and
        // an exact match passes even with a tolerance of zero.
        MALLOY_CHECK_EQ(quiet_near(1.0, 1.0 + 1e-12, 1e-9), 0);
        MALLOY_CHECK_EQ(quiet_near(1.0, 1.2, 1e-9), 1);
        MALLOY_CHECK_EQ(quiet_near(2.5, 2.5, 0.0), 0);

        // A NaN must FAIL rather than pass. Every comparison involving a NaN is
        // false, so the natural form (difference > epsilon) is false for a NaN
        // and silently accepts it. Two of these use an enormous tolerance to
        // show the rejection is about the NaN and not about the magnitude.
        MALLOY_CHECK_EQ(quiet_near(nan, 0.0, 0.0), 1);
        MALLOY_CHECK_EQ(quiet_near(0.0, nan, 1e9), 1);
        MALLOY_CHECK_EQ(quiet_near(nan, nan, 1e9), 1);

        // A NaN tolerance is rejected too, so a bad epsilon cannot quietly
        // accept every value handed to it.
        MALLOY_CHECK_EQ(quiet_near(1.0, 1.0, nan), 1);

        // Infinity minus infinity is NaN, so two infinities are not near each
        // other. Use MALLOY_CHECK_EQ for a value expected to be infinite.
        MALLOY_CHECK_EQ(quiet_near(inf, inf, 0.0), 1);
        MALLOY_CHECK_EQ(quiet_near(inf, 1.0, 1e9), 1);
    }

    // --- The same per component for the Vec2 form: a NaN in either component
    //     on its own is enough to fail. ---
    {
        MALLOY_CHECK_EQ(quiet_vec2_near(XY{1.0, 2.0}, XY{1.0, 2.0}, 0.0), 0);
        MALLOY_CHECK_EQ(quiet_vec2_near(XY{1.0, 2.0}, XY{1.0, 2.5}, 1e-9), 1);
        MALLOY_CHECK_EQ(quiet_vec2_near(XY{nan, 2.0}, XY{0.0, 2.0}, 1e9), 1);
        MALLOY_CHECK_EQ(quiet_vec2_near(XY{1.0, nan}, XY{1.0, 0.0}, 1e9), 1);
    }

    // --- A failure has to be REPORTED, not merely counted. A check that failed
    //     for some unrelated reason would otherwise look identical. ---
    {
        std::string reported;
        MALLOY_CHECK_EQ(quiet_near(nan, 0.0, 1.0, &reported), 1);
        MALLOY_CHECK_TRUE(reported.find("MALLOY_CHECK_NEAR failed") !=
                          std::string::npos);
    }

    // --- And a passing check says nothing at all. ---
    {
        std::string reported;
        MALLOY_CHECK_EQ(quiet_near(3.0, 3.0, 1e-9, &reported), 0);
        MALLOY_CHECK_TRUE(reported.empty());
    }

    std::cout << "malloy_smoke_tests passed\n";
    return 0;
}
