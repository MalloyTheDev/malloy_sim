#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <malloy/view3d/view3d.hpp>

#include "picture.hpp"
#include "window.hpp"

// The pixel self-test (`--self-test`): the window draws the frame's subjects
// flat on an exact background with no multisampling, reads the back buffer, and
// this checks the picture against the render model's own projection. A
// convention slip between the model and the drawing (a mirrored axis, rows read
// bottom-up, radians handed to a degrees API, the aspect inverted, y up) moves
// the picture away from where the model says it is, and fails.
namespace viewer
{
struct PixelReport
{
    std::size_t drawn_pixels{0}; // pixels that are not the background
    std::size_t subject_probes{0};
    std::size_t clear_probes{0};
    std::string failure; // empty when every check passed

    bool ok() const { return failure.empty(); }
};

// (a) the picture is not blank; (b) the pixel at the model's projection of
// every subject at least 2 px across is drawn; (c) every pixel of a 16 px grid
// whose view ray clears every subject's bounding sphere by more than
// 3 px + 0.1 of its radius is exactly the background.
PixelReport check_pixels(const std::vector<std::uint8_t>& rgba, int width, int height,
                         Rgb background, const malloy::view3d::Frame& frame, const Look& look);
} // namespace viewer
