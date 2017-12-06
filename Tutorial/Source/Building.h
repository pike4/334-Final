#pragma once
#include <string>


class Building 
{
public:
	int numSpires;
	double height;
	double width;
	double length;


	static Building prototype;
	static void genPrototype(std::string rules);
};