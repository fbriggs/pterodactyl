// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include "pImageBytes.h"
#include <vector>
using namespace std;

class pDrawPrimitives
{
public:
	
	static void ellipse_filled_1Channel(pImageBytes& img, float cx, float cy, float rA, float rB, float angle, int color)
	{
		for(float alpha = 0; alpha <= 1.0; alpha += 0.001)
			ellipse_1Channel(img, cx, cy, alpha * rA, alpha * rB, angle, color);
	}
	
	static void ellipse_1Channel(pImageBytes& img, float cx, float cy, float rA, float rB, float angle, int color)
	{
		vector< pair<float, float> > points;
		
		for(float t = 0; t < 2 * M_PI; t += 2 * M_PI / 100 )
		{
			float x = cx + rA * cos(t) * cos(angle) - rB * sin(t) * sin(angle);
			float y = cy + rA * cos(t) * sin(angle) + rB * sin(t) * cos(angle);
			
			points.push_back(pair<float,float>(x,y));
		}
		
		for(int i = 1; i < points.size(); ++i)
			line_1Channel(img, points[i].first, points[i].second, points[i-1].first, points[i-1].second, color);
		
		line_1Channel(img, points[0].first, points[0].second, points[points.size() - 1].first, points[points.size() - 1].second, color);
	}
	
	static void ellipse(pImageBytes& img, float cx, float cy, float rA, float rB, float angle, pRGB color)
	{
		vector< pair<float, float> > points;

		for(float t = 0; t < 2 * M_PI; t += 2 * M_PI / 100 )
		{
			float x = cx + rA * cos(t) * cos(angle) - rB * sin(t) * sin(angle);
			float y = cy + rA * cos(t) * sin(angle) + rB * sin(t) * cos(angle);
			
			points.push_back(pair<float,float>(x,y));
		}
		
		for(int i = 1; i < points.size(); ++i)
			line(img, points[i].first, points[i].second, points[i-1].first, points[i-1].second, color);
		
		line(img, points[0].first, points[0].second, points[points.size() - 1].first, points[points.size() - 1].second, color);
	}
	
	static void boxOutline(pImageBytes& img, int minX, int maxX, int minY, int maxY, float val, float alpha, int channel = 0)
	{
		for(int i = minX; i <= maxX; ++i)
		{
			img.setPixel(i, minY, channel, img.getPixel(i, minY, channel) * (1-alpha) + alpha * val);
			img.setPixel(i, maxY, channel, img.getPixel(i, maxY, channel) * (1-alpha) + alpha * val);
		}
		
		for(int j = minY +1; j < maxY; ++j)
		{
			img.setPixel(minX, j, channel, img.getPixel(minX, j, channel) * (1-alpha) + alpha * val);
			img.setPixel(maxX, j, channel, img.getPixel(maxX, j, channel) * (1-alpha) + alpha * val);
		}
	}
	
	static void boxOutline(pImageBytes& img, int minX, int maxX, int minY, int maxY, pRGB color)
	{
		for(int i = minX; i <= maxX; ++i)
		{
			img.setPixelRGB(i, minY, color);
			img.setPixelRGB(i, maxY, color);
		}
		
		for(int j = minY +1; j < maxY; ++j)
		{
			img.setPixelRGB(minX, j, color);
			img.setPixelRGB(maxX, j, color);
		}
	}
	
	static void circle(pImageBytes& img, float xc, float yc, float r, pRGB color)
	{
		float delta = 0.05;
		for(float theta = 0; theta <= 2 * M_PI; theta += delta)
		{
			float x1 = xc + r * cos(theta);
			float y1 = yc + r * sin(theta);
			
			float x2 = xc + r * cos(theta+delta);
			float y2 = yc + r * sin(theta+delta);
			
			line(img, x1, y1, x2, y2, color);
		}
	}
	
	
	static void line(pImageBytes& img, int x1, int y1,  int x2, int y2, pRGB color)
	{
		int delta_x = std::abs(x2 - x1) << 1;
		int delta_y = std::abs(y2 - y1) << 1;
		
		// if x1 == x2 or y1 == y2, then it does not matter what we set here
		signed char ix = x2 > x1?1:-1;
		signed char iy = y2 > y1?1:-1;
		
		img.setPixelRGB_(x1, y1, color);
		
		if (delta_x >= delta_y)
		{
			// error may go below zero
			int error = delta_y - (delta_x >> 1);
			
			while (x1 != x2)
			{
				if (error >= 0)
				{
					if (error || (ix > 0))
					{
						y1 += iy;
						error -= delta_x;
					}
					// else do nothing
				}
				// else do nothing
				
				x1 += ix;
				error += delta_y;
				
				img.setPixelRGB_(x1, y1, color);
			}
		}
		else
		{
			// error may go below zero
			int error = delta_x - (delta_y >> 1);
			
			while (y1 != y2)
			{
				if (error >= 0)
				{
					if (error || (iy > 0))
					{
						x1 += ix;
						error -= delta_y;
					}
					// else do nothing
				}
				// else do nothing
				
				y1 += iy;
				error += delta_x;
				
				img.setPixelRGB_(x1, y1, color);
			}
		}
	}
	
	
	static void line_1Channel(pImageBytes& img, int x1, int y1,  int x2, int y2, int color)
	{
		int delta_x = std::abs(x2 - x1) << 1;
		int delta_y = std::abs(y2 - y1) << 1;
		
		// if x1 == x2 or y1 == y2, then it does not matter what we set here
		signed char ix = x2 > x1?1:-1;
		signed char iy = y2 > y1?1:-1;
		
		img.setPixel_(x1, y1, 0, color);
		
		if (delta_x >= delta_y)
		{
			// error may go below zero
			int error = delta_y - (delta_x >> 1);
			
			while (x1 != x2)
			{
				if (error >= 0)
				{
					if (error || (ix > 0))
					{
						y1 += iy;
						error -= delta_x;
					}
					// else do nothing
				}
				// else do nothing
				
				x1 += ix;
				error += delta_y;
				
				img.setPixel_(x1, y1, 0, color);
			}
		}
		else
		{
			// error may go below zero
			int error = delta_x - (delta_y >> 1);
			
			while (y1 != y2)
			{
				if (error >= 0)
				{
					if (error || (iy > 0))
					{
						x1 += ix;
						error -= delta_y;
					}
					// else do nothing
				}
				// else do nothing
				
				y1 += iy;
				error += delta_x;
				
				img.setPixel_(x1, y1, 0, color);
			}
		}
	}
	
	
};
