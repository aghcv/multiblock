#include "fastvessels/io_utils.hpp"
#include <cassert>

int main() {
    try {
        auto poly = fastvessels::ReadPolyData("data/example.vtp");
        assert(poly->GetNumberOfPoints() > 0);
    } catch (...) {
        return 1;
    }
    return 0;
}
