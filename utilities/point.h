//
// Created by rosa on 11/19/24.
//

#ifndef POINT_H
#define POINT_H
#include <vector>  // For std::vector
#include <stdexcept>  // For std::invalid_argument


template<typename num_type>
class Point {
public:
    num_type x = 0;
    num_type y = 0;

    Point() = default;

    Point(num_type, num_type);

    bool operator<=(const Point &other) const {
        return (x <= other.x) && (y <= other.y);
    }
    Point operator*(double c) const {
        return Point(x*c, y*c);
    }
    Point operator+(const Point& other) const {
        return Point{x + other.x, y + other.y};
    }
    Point operator-(const Point& other) const {
        return Point{x - other.x, y - other.y};
    }
    bool operator != (const Point& other) const {
        return (x != other.x || y != other.y);
    }
    friend Point operator*(double scalar, const Point& point) {
        return Point(point.x * scalar, point.y * scalar);
    }
    double norm() {
        return std::sqrt(x*x + y*y);
    }
};


template<typename num_type>
Point<num_type>::Point(const num_type xCoord, const num_type yCoord)
    : x(xCoord), y(yCoord) {
}

/**
 * Creates a point with the minimum x and y coordinates from a vector of points.
 *
 * @tparam num_type The numeric type used for the coordinates.
 * @param points A vector of Point objects.
 * @return The Point object with the minimum x and y coordinates.
 * @throws std::invalid_argument if the points vector is empty.
 */
template<typename num_type>
Point<num_type> getMinPoint(const std::vector<Point<num_type> > &points) {
    if (points.empty()) {
        throw std::invalid_argument("The points vector is empty");
    }

    Point<num_type> minPoint = points[0];

    for (const auto &point: points) {
        if (point.x < minPoint.x) {
            minPoint.x = point.x;
        }
        if (point.y < minPoint.y) {
            minPoint.y = point.y;
        }
    }

    return minPoint;
}

/**
 * Creates a point with the maximum x and y coordinates from a vector of points.
 *
 * @tparam num_type The numeric type used for the coordinates.
 * @param points A vector of Point objects.
 * @return The Point object with the maximum x and y coordinates.
 * @throws std::invalid_argument if the points vector is empty.
 */

template<typename num_type>
Point<num_type> getMaxPoint(const std::vector<Point<num_type> > &points) {
    if (points.empty()) {
        throw std::invalid_argument("The points vector is empty");
    }

    Point<num_type> maxPoint = points[0];

    for (const auto &point: points) {
        if (point.x > maxPoint.x) {
            maxPoint.x = point.x;
        }
        if (point.y > maxPoint.y) {
            maxPoint.y = point.y;
        }
    }

    return maxPoint;
}

#endif //POINT_H

