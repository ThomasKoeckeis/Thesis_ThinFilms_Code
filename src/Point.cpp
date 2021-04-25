/*
 * Point.cpp
 *
 *  Created on: 25.03.2020
 *      Author: Köcki
 */

#include "Point.h"

Point::Point() {
	// TODO Auto-generated constructor stub
	x = 0;
	y = 0;
	z = 0;

}

Point::Point(float X,float Y,float Z) {
	// TODO Auto-generated constructor stub
	x = X;
	y = Y;
	z = Z;

}

Point::~Point() {
	// TODO Auto-generated destructor stub
}

Point::Point(const Point &other) {
	// TODO Auto-generated constructor stub
	x = other.x;
	y = other.y;
	z = other.z;
}

