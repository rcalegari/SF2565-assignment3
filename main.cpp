#include <iostream>
#include "utilities/grid.h"
#include <chrono>
#include <fstream>

void testPerformance(const int N, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::out | std::ios::trunc);
    for (int n = 3; n <= N; n++ ) {
        auto start = std::chrono::high_resolution_clock::now();

        Domain domain;
        Grid grid(n);
        grid.createGrid(domain);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time = end - start;
        std::cout << "Elapsed time for n = "<<n<<" : " << time.count() << " seconds\n";
        outFile << n << " " << time.count() << "\n"; // Write n and time to file
    }
    outFile.close();
}


int main() {
    // check if BottomCurve implementation works
    BottomCurve bottom;
    const std::string filename = "bottom_curve_coordinates.txt";
    bottom.writeCoordinatesToFile(filename);

    // DEFAULT DOMAIN //
    // the domain is defined as requested in the assignment
    Domain domain;
    int n = 50;
    Grid grid(n);
    // grid asked for in the assignment is saved in xGrid.txt and yGrid.txt
    grid.createGrid(domain);
    grid.printToFile("xGrid.txt", "yGrid.txt");

    // USER-DEFINED DOMAIN //
    // The domain is made up by the left, right and top boundaries being straight lines
    // and the bottom line being a smooth curve. The user has the freedom to choose the corners
    // of the domain and the equation of the curve.
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
    // grid defined by user is saved in xGrid_new.txt and yGrid_new.txt
    grid2.createGrid(domain2);
    grid2.printToFile("xGrid_new.txt", "yGrid_new.txt");

    // COMMENT OUT TO WORK WITH CACHE
    testPerformance(40, "performanceWithoutCache");
    // COMMENT OUT TO WORK WITHOUT CACHE
    // testPerformance(40, "performanceCache");

    return 0;
}

