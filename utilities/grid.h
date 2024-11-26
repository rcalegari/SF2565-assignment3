//
// Created by rosa on 11/22/24.
//

#ifndef GRID_H
#define GRID_H
#include<Eigen/Eigen>
#include<iostream>
#include<fstream>
#include "curve.h"
#include "point.h"

class Grid {
public:
    Eigen::MatrixXd xGridPoints;
    Eigen::MatrixXd yGridPoints;
    int n;

    Grid(int const n)
        : n(n), xGridPoints(n, n), yGridPoints(n, n) {
        if (n < 3) {
            throw std::invalid_argument("Error: n must be bigger than 2.");
        }
    };

    void createGrid( Domain& domain, const std::string& filenameX, const std::string& filenameY) {
        double hx = (domain.bottomBoundary.at(1).x - domain.bottomBoundary.at(0).x) / n;
        double hy = (domain.leftBoundary.at(1).y - domain.leftBoundary.at(0).y) / n;
        double h = 1./(n - 1);
        for (int j = 0; j < n; j++) {
            double eta = j * h;
            for (int i = 0; i < n; i++) {
                double xi = i * h;

                Point<double> bottomPoint = domain.bottomBoundary.at(xi);
                Point<double> topPoint = domain.topBoundary.at(xi);
                Point<double> leftPoint = domain.leftBoundary.at(eta);
                Point<double> rightPoint = domain.rightBoundary.at(eta);
                // std::cout << "  Bottom point: (" << bottomPoint.x << ", " << bottomPoint.y << ")" << std::endl;
                // std::cout << "  Top point: (" << topPoint.x << ", " << topPoint.y << ")" << std::endl;
                // std::cout << "  Left point: (" << leftPoint.x << ", " << leftPoint.y << ")" << std::endl;
                // std::cout<< "Right point: (" << rightPoint.x << ", " << rightPoint.y << ")" << std::endl;

                Point<double> gridPoint = 0.5 * (
                    (1 - eta) * bottomPoint + eta * topPoint +
                    (1 - xi) * leftPoint + xi * rightPoint ) ; //-( (1 - eta) * leftPoint + eta * rightPoint);
                xGridPoints(i, j) = gridPoint.x;
                yGridPoints(i, n-1-j) = gridPoint.y;

                // std::cout<< "grid point: (" << gridPoint.x << ", " << gridPoint.y << ")" << std::endl;
            }
        }
        std::ofstream xFile(filenameX, std::ios::out | std::ios::trunc);
        std::ofstream yFile(filenameY, std::ios::out | std::ios::trunc);
        if (!xFile.is_open() || !yFile.is_open()) {
            throw std::ios_base::failure("Failed to open output files");
        }
        // for (int row =0; row <= n; row++) {
        //     for (int col = 0; col <= n; col++) {
        for (int col = 0; col < n; col++) {
            for (int row =0; row < n; row++) {
                xFile<<xGridPoints(row, col)<<" ";
                yFile<<yGridPoints(row, col)<<" ";
            }
            xFile<<"\n";
            yFile<<"\n";
        }
        xFile.close();
        yFile.close();
    }
};

#endif //GRID_H