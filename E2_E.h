#pragma once
#include <math.h>
#include <stdint.h>
#include "gcclib.h"

static const uint32_t numCoeffs = 2;

template<typename T> struct GCCLIB_API E2_EMV
{
    T coeffs[numCoeffs];

    E2_EMV<T>() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }

    E2_EMV<T>(T S, T E0E1)
    {
        coeffs[0] = S;
        coeffs[1] = E0E1;
    }

    E2_EMV<T>& operator=(const E2_EMV<T>& rhs)
    {
        if (this != &rhs)
            memcpy(&coeffs[0], &rhs.coeffs[0], numCoeffs * sizeof(T));

        return *this;
    }

    E2_EMV<T>& operator+=(const E2_EMV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] += rhs.coeffs[i];

        return *this;
    }

    E2_EMV<T>& operator-=(const E2_EMV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] -= rhs.coeffs[i];

        return *this;
    }

    E2_EMV<T>& operator|=(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] - tmp.coeffs[1] * rhs.coeffs[1];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0];

        return *this;
    }

    E2_EMV<T>& operator^=(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;


        return *this;
    }

    E2_EMV<T>& operator*=(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] - tmp.coeffs[1] * rhs.coeffs[1];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0];

        return *this;
    }

    E2_EMV<T>& operator*=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] *= rhs;

        return *this;
    }

    E2_EMV<T>& operator/=(const E2_EMV<T>& rhs)
    {
        E2_EMV<T>tmp = ~rhs;
        T mag = rhs.Mag();
        (*this) *= tmp;
        (*this) /= (mag * mag);
        return *this;
    }

    E2_EMV<T>& operator/=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] /= rhs;

        return *this;
    }

    E2_EMV<T> operator+(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E2_EMV<T> operator-(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp -= rhs;
        return tmp;
    }

    E2_EMV<T> operator|(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp |= rhs;
        return tmp;
    }

    E2_EMV<T> operator^(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp ^= rhs;
        return tmp;
    }

    E2_EMV<T> operator*(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp *= rhs;
        return tmp;
    }

    E2_EMV<T> operator*(const T rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E2_EMV<T> operator/(const E2_EMV<T>& rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E2_EMV<T> operator/(const T rhs)
    {
        E2_EMV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E2_EMV<T> operator~()
    {
        E2_EMV<T> tmp = *this;
        tmp.coeffs[1] = -tmp.coeffs[1];

        return tmp;
    }

    T S(){ return coeffs[0]; }
    void S(T val){ coeffs[0] = val; }

    T E0E1(){ return coeffs[1]; }
    void E0E1(T val){ coeffs[1] = val; }

    T Mag()
    {
        T mag = 0.0;
        for (uint32_t i = 0; i < nCoeffs; i++) mag += coeffs[i] * coeffs[i];
        return sqrt(mag);
    }

};