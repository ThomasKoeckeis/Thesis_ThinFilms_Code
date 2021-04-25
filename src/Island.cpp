/*
 * Island.cpp
 *
 *  Created on: 25.03.2020
 *      Author: Köcki
 */

#include "Island.h"
#include "Polygon.h"
#include "Point.h"
#include <list>

Island::Island() {
	// TODO Auto-generated constructor stub
	std::vector<Polygon> polygons = std::vector<Polygon>();
}

Island::~Island() {
	// TODO Auto-generated destructor stub
}

Island::Island(const Island &other) {
	// TODO Auto-generated constructor stub
	//for ( int i = 0; i < polygons.size; i++)
	//{
	//	polygons.push_back(other.polygons[i]);
	//}
		polygons.push_back(other.polygons[0]);
}

