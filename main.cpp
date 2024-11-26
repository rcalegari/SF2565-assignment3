#include <iostream>
// #include "utilities/curve.h"
#include "utilities/grid.h"
// #include "utilities/point.h"

int main() {
    Point<double> p0 = {-10, 0};
    Point<double> p1 = {5, 0};

    StraightLine line(p0, p1);
    line.print();

    // Check if at function works
    Point<double> p_mid = line.at(0.5);  // Midpoint
    std::cout << "Midpoint: (" << p_mid.x << ", " << p_mid.y << ")" << std::endl;
    Point<double> p_start = line.at(0.0);  // Start point
    std::cout << "Start point: (" << p_start.x << ", " << p_start.y << ")" << std::endl;
    Point<double> p_end = line.at(1.0);  // End point
    std::cout << "End point: (" << p_end.x << ", " << p_end.y << ")" << std::endl;

    // check if BottomCurve implementation works
    BottomCurve bottom;
    const std::string filename = "bottom_curve_coordinates.txt";
    bottom.writeCoordinatesToFile(filename);

    // check for default constructors
    Domain domain;
    int n = 10;
    Grid grid(n);
    grid.createGrid(domain, "xGrid.txt", "yGrid.txt");


    return 0;
}

