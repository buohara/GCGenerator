#pragma once
#include <math.h>
#include <stdint.h>
#include "gcclib.h"

static const uint32_t numCoeffs = 8;

template<typename T> struct GCCLIB_API E3MV
{
    T coeffs[numCoeffs];

    E3MV<T>() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }

    E3MV<T>(T S, T E0, T E1, T E2, T E0E1, T E0E2, T E1E2, T E0E1E2)
    {
        coeffs[0] = S;
        coeffs[1] = E0;
        coeffs[2] = E1;
        coeffs[3] = E2;
        coeffs[4] = E0E1;
        coeffs[5] = E0E2;
        coeffs[6] = E1E2;
        coeffs[7] = E0E1E2;
    }

    E3MV<T>& operator=(const E3MV<T>& rhs)
    {
        if (this != &rhs)
            memcpy(&coeffs[0], &rhs.coeffs[0], numCoeffs * sizeof(T));

        return *this;
    }

    E3MV<T>& operator+=(const E3MV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] += rhs.coeffs[i];

        return *this;
    }

    E3MV<T>& operator-=(const E3MV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] -= rhs.coeffs[i];

        return *this;
    }

    E3MV<T>& operator|=(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] + tmp.coeffs[1] * rhs.coeffs[1] + tmp.coeffs[2] * rhs.coeffs[2] + tmp.coeffs[3] * rhs.coeffs[3] - tmp.coeffs[4] * rhs.coeffs[4] - tmp.coeffs[5] * rhs.coeffs[5] - tmp.coeffs[6] * rhs.coeffs[6] - tmp.coeffs[7] * rhs.coeffs[7];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0] - tmp.coeffs[2] * rhs.coeffs[4] - tmp.coeffs[3] * rhs.coeffs[5] + tmp.coeffs[4] * rhs.coeffs[2] + tmp.coeffs[5] * rhs.coeffs[3] - tmp.coeffs[6] * rhs.coeffs[7] - tmp.coeffs[7] * rhs.coeffs[6];
        coeffs[2] = tmp.coeffs[0] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[4] + tmp.coeffs[2] * rhs.coeffs[0] - tmp.coeffs[3] * rhs.coeffs[6] - tmp.coeffs[4] * rhs.coeffs[1] + tmp.coeffs[5] * rhs.coeffs[7] + tmp.coeffs[6] * rhs.coeffs[3] + tmp.coeffs[7] * rhs.coeffs[5];
        coeffs[3] = tmp.coeffs[0] * rhs.coeffs[3] + tmp.coeffs[1] * rhs.coeffs[5] + tmp.coeffs[2] * rhs.coeffs[6] + tmp.coeffs[3] * rhs.coeffs[0] - tmp.coeffs[4] * rhs.coeffs[7] - tmp.coeffs[5] * rhs.coeffs[1] - tmp.coeffs[6] * rhs.coeffs[2] - tmp.coeffs[7] * rhs.coeffs[4];
        coeffs[4] = tmp.coeffs[0] * rhs.coeffs[4] + tmp.coeffs[3] * rhs.coeffs[7] + tmp.coeffs[4] * rhs.coeffs[0] - tmp.coeffs[5] * rhs.coeffs[6] + tmp.coeffs[6] * rhs.coeffs[5] + tmp.coeffs[7] * rhs.coeffs[3];
        coeffs[5] = tmp.coeffs[0] * rhs.coeffs[5] - tmp.coeffs[2] * rhs.coeffs[7] + tmp.coeffs[4] * rhs.coeffs[6] + tmp.coeffs[5] * rhs.coeffs[0] - tmp.coeffs[6] * rhs.coeffs[4] - tmp.coeffs[7] * rhs.coeffs[2];
        coeffs[6] = tmp.coeffs[0] * rhs.coeffs[6] + tmp.coeffs[1] * rhs.coeffs[7] - tmp.coeffs[4] * rhs.coeffs[5] + tmp.coeffs[5] * rhs.coeffs[4] + tmp.coeffs[6] * rhs.coeffs[0] + tmp.coeffs[7] * rhs.coeffs[1];
        coeffs[7] = tmp.coeffs[0] * rhs.coeffs[7] + tmp.coeffs[7] * rhs.coeffs[0];

        return *this;
    }

    E3MV<T>& operator^=(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;

        coeffs[4] = tmp.coeffs[1] * rhs.coeffs[2] - tmp.coeffs[2] * rhs.coeffs[1];
        coeffs[5] = tmp.coeffs[1] * rhs.coeffs[3] - tmp.coeffs[3] * rhs.coeffs[1];
        coeffs[6] = tmp.coeffs[2] * rhs.coeffs[3] - tmp.coeffs[3] * rhs.coeffs[2];
        coeffs[7] = tmp.coeffs[1] * rhs.coeffs[6] - tmp.coeffs[2] * rhs.coeffs[5] + tmp.coeffs[3] * rhs.coeffs[4] + tmp.coeffs[4] * rhs.coeffs[3] - tmp.coeffs[5] * rhs.coeffs[2] + tmp.coeffs[6] * rhs.coeffs[1];

        return *this;
    }

    E3MV<T>& operator*=(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] + tmp.coeffs[1] * rhs.coeffs[1] + tmp.coeffs[2] * rhs.coeffs[2] + tmp.coeffs[3] * rhs.coeffs[3] - tmp.coeffs[4] * rhs.coeffs[4] - tmp.coeffs[5] * rhs.coeffs[5] - tmp.coeffs[6] * rhs.coeffs[6] - tmp.coeffs[7] * rhs.coeffs[7];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0] - tmp.coeffs[2] * rhs.coeffs[4] - tmp.coeffs[3] * rhs.coeffs[5] + tmp.coeffs[4] * rhs.coeffs[2] + tmp.coeffs[5] * rhs.coeffs[3] - tmp.coeffs[6] * rhs.coeffs[7] - tmp.coeffs[7] * rhs.coeffs[6];
        coeffs[2] = tmp.coeffs[0] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[4] + tmp.coeffs[2] * rhs.coeffs[0] - tmp.coeffs[3] * rhs.coeffs[6] - tmp.coeffs[4] * rhs.coeffs[1] + tmp.coeffs[5] * rhs.coeffs[7] + tmp.coeffs[6] * rhs.coeffs[3] + tmp.coeffs[7] * rhs.coeffs[5];
        coeffs[3] = tmp.coeffs[0] * rhs.coeffs[3] + tmp.coeffs[1] * rhs.coeffs[5] + tmp.coeffs[2] * rhs.coeffs[6] + tmp.coeffs[3] * rhs.coeffs[0] - tmp.coeffs[4] * rhs.coeffs[7] - tmp.coeffs[5] * rhs.coeffs[1] - tmp.coeffs[6] * rhs.coeffs[2] - tmp.coeffs[7] * rhs.coeffs[4];
        coeffs[4] = tmp.coeffs[0] * rhs.coeffs[4] + tmp.coeffs[3] * rhs.coeffs[7] + tmp.coeffs[4] * rhs.coeffs[0] - tmp.coeffs[5] * rhs.coeffs[6] + tmp.coeffs[6] * rhs.coeffs[5] + tmp.coeffs[7] * rhs.coeffs[3] + tmp.coeffs[1] * rhs.coeffs[2] - tmp.coeffs[2] * rhs.coeffs[1];
        coeffs[5] = tmp.coeffs[0] * rhs.coeffs[5] - tmp.coeffs[2] * rhs.coeffs[7] + tmp.coeffs[4] * rhs.coeffs[6] + tmp.coeffs[5] * rhs.coeffs[0] - tmp.coeffs[6] * rhs.coeffs[4] - tmp.coeffs[7] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[3] - tmp.coeffs[3] * rhs.coeffs[1];
        coeffs[6] = tmp.coeffs[0] * rhs.coeffs[6] + tmp.coeffs[1] * rhs.coeffs[7] - tmp.coeffs[4] * rhs.coeffs[5] + tmp.coeffs[5] * rhs.coeffs[4] + tmp.coeffs[6] * rhs.coeffs[0] + tmp.coeffs[7] * rhs.coeffs[1] + tmp.coeffs[2] * rhs.coeffs[3] - tmp.coeffs[3] * rhs.coeffs[2];
        coeffs[7] = tmp.coeffs[0] * rhs.coeffs[7] + tmp.coeffs[7] * rhs.coeffs[0] + tmp.coeffs[1] * rhs.coeffs[6] - tmp.coeffs[2] * rhs.coeffs[5] + tmp.coeffs[3] * rhs.coeffs[4] + tmp.coeffs[4] * rhs.coeffs[3] - tmp.coeffs[5] * rhs.coeffs[2] + tmp.coeffs[6] * rhs.coeffs[1];

        return *this;
    }

    E3MV<T>& operator*=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] *= rhs;

        return *this;
    }

    E3MV<T>& operator/=(const E3MV<T>& rhs)
    {
        E3MV<T>tmp = ~rhs;
        T mag = rhs.Mag();
        (*this) *= tmp;
        (*this) /= (mag * mag);
        return *this;
    }

    E3MV<T>& operator/=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] /= rhs;

        return *this;
    }

    E3MV<T> operator+(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E3MV<T> operator-(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;
        tmp -= rhs;
        return tmp;
    }

    E3MV<T> operator|(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;
        tmp |= rhs;
        return tmp;
    }

    E3MV<T> operator^(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;
        tmp ^= rhs;
        return tmp;
    }

    E3MV<T> operator*(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;
        tmp *= rhs;
        return tmp;
    }

    E3MV<T> operator*(const T rhs)
    {
        E3MV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E3MV<T> operator/(const E3MV<T>& rhs)
    {
        E3MV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E3MV<T> operator/(const T rhs)
    {
        E3MV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E3MV<T> operator~()
    {
        E3MV<T> tmp = *this;
        tmp.coeffs[4] = -tmp.coeffs[4];
        tmp.coeffs[5] = -tmp.coeffs[5];
        tmp.coeffs[6] = -tmp.coeffs[6];
        tmp.coeffs[7] = -tmp.coeffs[7];

        return tmp;
    }

    T S(){ return coeffs[0]; }
    void S(T val){ coeffs[0] = val; }

    T E0(){ return coeffs[1]; }
    void E0(T val){ coeffs[1] = val; }

    T E1(){ return coeffs[2]; }
    void E1(T val){ coeffs[2] = val; }

    T E2(){ return coeffs[3]; }
    void E2(T val){ coeffs[3] = val; }

    T E0E1(){ return coeffs[4]; }
    void E0E1(T val){ coeffs[4] = val; }

    T E0E2(){ return coeffs[5]; }
    void E0E2(T val){ coeffs[5] = val; }

    T E1E2(){ return coeffs[6]; }
    void E1E2(T val){ coeffs[6] = val; }

    T E0E1E2(){ return coeffs[7]; }
    void E0E1E2(T val){ coeffs[7] = val; }

    T Mag()
    {
        T mag = 0.0;
        for (uint32_t i = 0; i < nCoeffs; i++) mag += coeffs[i] * coeffs[i];
        return sqrt(mag);
    }

};