#include <vector>
#include <map>
#include <stdio>
using namespace std;

enum STABasis
{
	e0,
	e1,
	e2,
	e3,
	e0e1,
	e0e2,
	e0e3,
	e1e2,
	e1e3,
	e2e3,
	e0e1e2,
	e0e1e3,
	e0e2e3,
	e1e2e3,
	e0e1e2e3,
};

struct STAMV
{
	float e0;
	float e1;
	float e2;
	float e3;
	float e0e1;
	float e0e2;
	float e0e3;
	float e1e2;
	float e1e3;
	float e2e3;
	float e0e1e2;
	float e0e1e3;
	float e0e2e3;
	float e1e2e3;
	float e0e1e2e3;
};

struct STARotor
{
	float e0e1;
	float e0e2;
	float e0e3;
	float e1e2;
	float e1e3;
	float e2e3;
	float e0e1e2e3;
};