#pragma once
#include <math.h>
#include <stdint.h>
#include "gcclib.h"

static const uint32_t numCoeffs = 4;

template<typename T> struct GCCLIB_API E3_EMV
{
    T coeffs[numCoeffs];

    E3_EMV<T>() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }

    E3_EMV<T>(T S, T E0E1, T E0E2, T E1E2)
    {
        coeffs[0] = S;
        coeffs[1] = E0E1;
        coeffs[2] = E0E2;
        coeffs[3] = E1E2;
    }

    E3_EMV<T>& operator=(const E3_EMV<T>& rhs)
    {
        if (this != &rhs)
            memcpy(&coeffs[0], &rhs.coeffs[0], numCoeffs * sizeof(T));

        return *this;
    }

    E3_EMV<T>& operator+=(const E3_EMV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] += rhs.coeffs[i];

        return *this;
    }

    E3_EMV<T>& operator-=(const E3_EMV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] -= rhs.coeffs[i];

        return *this;
    }

    E3_EMV<T>& operator|=(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] - tmp.coeffs[1] * rhs.coeffs[1] - tmp.coeffs[2] * rhs.coeffs[2] - tmp.coeffs[3] * rhs.coeffs[3];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0] - tmp.coeffs[2] * rhs.coeffs[3] + tmp.coeffs[3] * rhs.coeffs[2];
        coeffs[2] = tmp.coeffs[0] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[3] + tmp.coeffs[2] * rhs.coeffs[0] - tmp.coeffs[3] * rhs.coeffs[1];
        coeffs[3] = tmp.coeffs[0] * rhs.coeffs[3] - tmp.coeffs[1] * rhs.coeffs[2] + tmp.coeffs[2] * rhs.coeffs[1] + tmp.coeffs[3] * rhs.coeffs[0];

        return *this;
    }

    E3_EMV<T>& operator^=(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;


        return *this;
    }

    E3_EMV<T>& operator*=(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] - tmp.coeffs[1] * rhs.coeffs[1] - tmp.coeffs[2] * rhs.coeffs[2] - tmp.coeffs[3] * rhs.coeffs[3];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0] - tmp.coeffs[2] * rhs.coeffs[3] + tmp.coeffs[3] * rhs.coeffs[2];
        coeffs[2] = tmp.coeffs[0] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[3] + tmp.coeffs[2] * rhs.coeffs[0] - tmp.coeffs[3] * rhs.coeffs[1];
        coeffs[3] = tmp.coeffs[0] * rhs.coeffs[3] - tmp.coeffs[1] * rhs.coeffs[2] + tmp.coeffs[2] * rhs.coeffs[1] + tmp.coeffs[3] * rhs.coeffs[0];

        return *this;
    }

    E3_EMV<T>& operator*=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] *= rhs;

        return *this;
    }

    E3_EMV<T>& operator/=(const E3_EMV<T>& rhs)
    {
        E3_EMV<T>tmp = ~rhs;
        T mag = rhs.Mag();
        (*this) *= tmp;
        (*this) /= (mag * mag);
        return *this;
    }

    E3_EMV<T>& operator/=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] /= rhs;

        return *this;
    }

    E3_EMV<T> operator+(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E3_EMV<T> operator-(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp -= rhs;
        return tmp;
    }

    E3_EMV<T> operator|(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp |= rhs;
        return tmp;
    }

    E3_EMV<T> operator^(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp ^= rhs;
        return tmp;
    }

    E3_EMV<T> operator*(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp *= rhs;
        return tmp;
    }

    E3_EMV<T> operator*(const T rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E3_EMV<T> operator/(const E3_EMV<T>& rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E3_EMV<T> operator/(const T rhs)
    {
        E3_EMV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E3_EMV<T> operator~()
    {
        E3_EMV<T> tmp = *this;
        tmp.coeffs[1] = -tmp.coeffs[1];
        tmp.coeffs[2] = -tmp.coeffs[2];
        tmp.coeffs[3] = -tmp.coeffs[3];

        return tmp;
    }

    T S(){ return coeffs[0]; }
    void S(T val){ coeffs[0] = val; }

    T E0E1(){ return coeffs[1]; }
    void E0E1(T val){ coeffs[1] = val; }

    T E0E2(){ return coeffs[2]; }
    void E0E2(T val){ coeffs[2] = val; }

    T E1E2(){ return coeffs[3]; }
    void E1E2(T val){ coeffs[3] = val; }

    T Mag()
    {
        T mag = 0.0;
        for (uint32_t i = 0; i < nCoeffs; i++) mag += coeffs[i] * coeffs[i];
        return sqrt(mag);
    }

};