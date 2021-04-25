/*
 * Polygon.h
 *
 *  Created on: 25.03.2020
 *      Author: Köcki
 */

#ifndef POLYGON_H_
#define POLYGON_H_
#include "Point.h"
#include <vector>

 /// <summary>
 ///Class for Polygons with points and values.
 /// </summary>
class Polygon {
public:
	Polygon();
	virtual ~Polygon();
	Polygon(const Polygon &other);
	Polygon(int pointindex[], int size,float _stress_scalar);
	Polygon(int pointindex[], int size);

	std::vector<int> pointsindex;
	float stress_scalar;
	Point stress_vector;
};

#endif /* POLYGON_H_ */
