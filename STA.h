#include <vector>
#include <map>
#include <stdio>
#include <stdint.h>

using namespace std;

enum STABasis
{
	s,
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
	const uint32_t numCoeffs = 16;
	float coeffs[numCoeffs];

	STAMV& operator=(const STAMV& rhs)
	{
		if (this != &rhs)
		{
			for (uint32_t i = 0; i < numCoeffs; i++)
			{
				coeffs[i] = rhs.coeffs[i];
			}
		}
		return *this;
	}

	STAMV& operator+=(const STAMV& rhs)
	{
		for (uint32_t i = 0; i < numCoeffs; i++)
		{
			coeffs[i] += rhs.coeffs[i];
		}
		return *this;
	}

	STAMV& operator-=(const STAMV& rhs)
	{
		for (uint32_t i = 0; i < numCoeffs; i++)
		{
			coeffs[i] -= rhs.coeffs[i];
		}
		return *this;
	}

};

struct STARotor
{
	const uint32_t numCoeffs = 8.0;
	float coeffs[numCoeffs];

};

