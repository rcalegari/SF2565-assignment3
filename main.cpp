#include <iostream>
// #include "utilities/curve.h"
#include "utilities/grid.h"
// #include "utilities/point.h"

int main() {

    // check if BottomCurve implementation works
    BottomCurve bottom;
    const std::string filename = "bottom_curve_coordinates.txt";
    bottom.writeCoordinatesToFile(filename);

    // check for default domain
    Domain domain;
    int n = 10;
    Grid grid(n);
    grid.createGrid(domain, "xGrid.txt", "yGrid.txt");

    // check for user-defined domain
    Point<double> pp1 = {-5, 7};
    Point<double> pp2 = {4, 5};
    Point<double> pp3 = {-3, 11};
    Point<double> pp4 = {2, 12};

    // straight line
    // auto f = [](double x) -> double {
    //     return 53 / 9. - 2 * x / 9.;
    // };

    // curve
    auto f = [](double x) {
        return 7 + (5 - 7) * std::pow(std::sin(M_PI * (x + 5) / 18), 2);
    };

    StraightLine left(pp1, pp3);
    StraightLine top(pp3, pp4);
    StraightLine right(pp2, pp4);
    BottomCurve bottom2(f, pp1, pp2);
    Domain domain2(left, right, top, bottom2);

    Grid grid2(n);
    grid2.createGrid(domain2, "xGrid_new.txt", "yGrid_new.txt");

    return 0;
}

