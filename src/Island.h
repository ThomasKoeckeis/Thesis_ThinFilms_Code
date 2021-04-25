/*
 * Island.h
 *
 *  Created on: 25.03.2020
 *      Author: Köcki
 */

#ifndef ISLAND_H_
#define ISLAND_H_
#include "Polygon.h"
#include "Point.h"
#include <list>

 /// <summary>
 ///Class for Islands (unused).
 /// </summary>
class Island {
public:
	Island();
	virtual ~Island();
	Island(const Island &other);
	std::vector<Polygon> polygons;
};

#endif /* ISLAND_H_ */
