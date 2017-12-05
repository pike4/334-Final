#pragma once
#include <string>


class Building 
{
public:
	int numSpires;
	float height;
	float width;
	float length;


	static Building prototype;
	static void genPrototype(std::string rules);
};