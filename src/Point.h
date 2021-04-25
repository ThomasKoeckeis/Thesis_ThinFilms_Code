/*
 * Point.h
 *
 *  Created on: 25.03.2020
 *      Author: K�cki
 */

#ifndef POINT_H_
#define POINT_H_


 /// <summary>
 ///Class for 3D-points.
 /// </summary>
class Point {
public:
	Point();
	Point(float X, float Y, float Z);
	virtual ~Point();
	Point(const Point &other);
	float x;
	float y;
	float z;
private:                             // privat
   
};

#endif /* POINT_H_ */
