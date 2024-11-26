//
// Created by rosa on 11/19/24.
//

#ifndef CURVE_H
#define CURVE_H

#include<iostream>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <Eigen/Eigen>
#include "point.h"
#include "../../SF2565-Assignments/include/algo.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include<cmath>

class Curve {
public:
    virtual ~Curve() = default;
    virtual Point<double> at(double t) const = 0;
};

class StraightLine : public Curve {

public:
    Point<double> a;
    Point<double> b;
    StraightLine(const Point<double>& p0, const Point<double>& p1)
        : a(p0), b(p1){};
    Point<double> at(double t) const override {
        return a * (1 - t) + b * t;
    }
    void print() {
        std::cout << "StraightLine from (" << a.x << ", " << a.y
                  << ") to (" << b.x << ", " << b.y << ")" << std::endl;
    }
};

class EquationCurve : public Curve {
public:
    virtual ~EquationCurve() = default;

    Point<double> at(double t) const override {
        double sTot = arcLength(1.);
        double s = t * sTot;
        double tPrime = findTgivenArcLength(s);
        return gamma(tPrime);
    };
private:
    virtual Point<double> gamma(double t) const = 0;
    virtual Point<double> gammaprime(double t) const = 0;

    double arcLength( double t) const {
        boost::math::quadrature::gauss_kronrod<double, 15> integrator;
        auto normGammaPrime = [this](double u) {
            Point<double> derivative = gammaprime(u);
            return std::sqrt(derivative.x * derivative.x + derivative.y * derivative.y);
        };
        return integrator.integrate(normGammaPrime, 0., t, 1e-6, 1e-10);
    };

    // root-finding algorithm
    double findTgivenArcLength(double s) const {
        double t0 = 0.;
        double t1 = 1.;
        double tStar;
        while (true) {
            tStar = (t0 + t1) / 2.;
            if (t1 - t0 < 1e-5) {
                return tStar;
            }
            if ( arcLength(tStar) < s ) {
                t0 = tStar;
            }
            else {
                t1 = tStar;
            }
        }
    }

};

class BottomCurve : public EquationCurve {
private:
    std::function<double(double)> f;
    Point<double> leftPoint;
    Point<double> rightPoint;
    // these points have to coincide with the bottom points of
    // leftBoundary and rightBoundary
public:
    // the default function is the one given in the assignment
    BottomCurve()
        : leftPoint(Point<double>(-10, 0)),
        rightPoint(Point<double>(5, 0)),
        f([](double x) -> double {
            if (x >= -10 && x <= -3) {
                return 0.5 * 1 / (1 + std::exp(-3 * (x + 6)));
            } else if (x > -3 && x <= 5) {
                return 0.5 * 1 / (1 + std::exp(3 * x));
            } else {
                throw std::out_of_range("x is out of the boundary.");
            }
        }) {}

    BottomCurve(const std::function<double(double)>& func, Point<double> left, Point<double> right)
    : f(func), leftPoint(left), rightPoint(right) {
        if (f(leftPoint.x) != leftPoint.y || f(rightPoint.x) != rightPoint.y) {
            throw std::invalid_argument("f doesn't pass through the given corner points");
        };
    }

    Point<double> gamma(double t) const override {
        double x = (1 - t) * leftPoint.x + t * rightPoint.x;
        return Point<double>(x, f(x));
    };

    Point<double> gammaprime(double t) const override {
        double dx_dt = - rightPoint.x + leftPoint.x;
        double x = (1 - t) * leftPoint.x + t * rightPoint.x;
        // central finite difference
        double dx = 1e-7;
        double df_dx = (f(x + dx) - f(x - dx)) / (2 * dx);
        double df_dt = df_dx * dx_dt;
        return Point<double>(dx_dt, df_dt);
    };

    void printCoordinates() const {
        int num_points = 100; // Number of points to sample
        for (int i = 0; i <= num_points; ++i) {
            double t = static_cast<double>(i) / num_points;
            Point<double> point = gamma(t);
            std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
        }
    }
    void writeCoordinatesToFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::ios_base::failure("Failed to open file for writing.");
        }
        int num_points = 100; // Number of points to sample
        for (int i = 0; i <= num_points; ++i) {
            double t = static_cast<double>(i) / num_points;
            Point<double> point = gamma(t);
            file << point.x << " " << point.y << std::endl;
        }
        file.close();
    }


};

class Domain {
public:
    StraightLine leftBoundary;
    StraightLine rightBoundary;
    StraightLine topBoundary;
    BottomCurve bottomBoundary;
    // StraightLine bottomBoundary;

    // Domain()
    //     : leftBoundary(Point<double>(-10, 0), Point<double>(-10, 3)),
    //       rightBoundary(Point<double>(5, 0), Point<double>(5, 3)),
    //       topBoundary(Point<double>(-10, 3), Point<double>(5, 3)),
    //       bottomBoundary(Point<double>(-10, 0), Point<double>(5, 0)) {}
    Domain()
            : leftBoundary(Point<double>(-10, 0), Point<double>(-10, 3)),
              rightBoundary(Point<double>(5, 0), Point<double>(5, 3)),
              topBoundary(Point<double>(-10, 3), Point<double>(5, 3)),
              bottomBoundary() {};

    Domain(StraightLine& left, StraightLine& right, StraightLine& top, BottomCurve& bottom)
        : leftBoundary(left), rightBoundary(right), topBoundary(top), bottomBoundary(bottom) {
        // check if domain is closed
        // double norm = (left.at(0.) - bottom.at(0.)).norm();
        // std::cout<<norm<<std::endl;
        // std::cout<<"("<<left.at(0.).x<<", "<<left.at(0.).y<<")"<<std::endl;
        // std::cout<<"("<<bottom.at(0.).x<<", "<<bottom.at(0.).y<<")"<<std::endl;
        // Point<double> diff = left.at(0.) - bottom.at(0.);
        // std::cout << "Difference: (" << diff.x << ", " << diff.y << ")" << std::endl;

        double tol = 1e-4;
        if ((left.at(0.) - bottom.at(0.)).norm() > tol) {
            throw std::invalid_argument("The domain is not closed - pp1 doesn't coincide.");
        }
        if ((left.at(1.) - top.at(0.)).norm() > tol) {
            throw std::invalid_argument("The domain is not closed - pp3 doesn't coincide.");
        }
        if ((right.at(0.) - bottom.at(1.)).norm() > tol) {
            throw std::invalid_argument("The domain is not closed - pp2 doesn't coincide.");
        }
        if ((right.at(1.) - top.at(1.)).norm() > tol) {
            throw std::invalid_argument("The domain is not closed - pp4 doesn't coincide.");
        }
    };

    Domain(const BottomCurve& curve)
        : leftBoundary(Point<double>(-10, 0), Point<double>(-10, 3)),
          rightBoundary(Point<double>(5, 0), Point<double>(5, 3)),
          topBoundary(Point<double>(-10, 3), Point<double>(5, 3)),
          bottomBoundary(curve) {}
};



#endif //CURVE_H
