/*
 * Polygon.cpp
 *
 *  Created on: 25.03.2020
 *      Author: Köcki
 */

#include "Polygon.h"
#include "Point.h"
#include <vector>

Polygon::Polygon() {
	// TODO Auto-generated constructor stub
	std::vector<int> pointsindex =  std::vector<int>();
	stress_scalar = 0;
	Point stress_vector = Point::Point();
}

Polygon::~Polygon() {
	// TODO Auto-generated destructor stub
}

Polygon::Polygon(const Polygon &other) {
	// TODO Auto-generated constructor stub
	pointsindex = other.pointsindex;
	stress_scalar = other.stress_scalar;
	stress_vector = other.stress_vector;
}

Polygon::Polygon(int pointindex[],int size, float _stress_scalar) {
	// TODO Auto-generated constructor stub
	for (int i = 0; i < size; i++) {
		pointsindex.push_back( pointindex[i]);
	}
	stress_scalar = _stress_scalar;
	stress_vector = Point(0,0,0);
}

Polygon::Polygon(int pointindex[], int size) {
	// TODO Auto-generated constructor stub
	for (int i = 0; i < size; i++) {
		pointsindex.push_back(pointindex[i]);
	}
	stress_scalar = 0;
	stress_vector = Point(0, 0, 0);
}
