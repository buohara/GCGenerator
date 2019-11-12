static const uint32_t numCoeffs = 4;

template<typename T> struct E2MV
{
    T coeffs[numCoeffs];

    E2<T>() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }

    E2<T>(T S, T E0, T E1, T E0E1)
    {
        coeffs[0] = S;
        coeffs[1] = E0;
        coeffs[2] = E1;
        coeffs[3] = E0E1;
    }

    E2MV<T>& operator=(const E2MV<T>& rhs)
    {
        if (this != &rhs)
            memcpy(&coeffs[0], &rhs.coeffs[0], numCoeffs * sizeof(T));

        return *this;
    }

    E2MV<T>& operator+=(const E2MV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] += rhs.coeffs[i];

        return *this;
    }

    E2MV<T>& operator-=(const E2MV<T>& rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] -= rhs.coeffs[i];

        return *this;
    }

    E2MV<T>& operator|=(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] + tmp.coeffs[1] * rhs.coeffs[1] + tmp.coeffs[2] * rhs.coeffs[2] - tmp.coeffs[3] * rhs.coeffs[3];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0] - tmp.coeffs[2] * rhs.coeffs[3] + tmp.coeffs[3] * rhs.coeffs[2];
        coeffs[2] = tmp.coeffs[0] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[3] + tmp.coeffs[2] * rhs.coeffs[0] - tmp.coeffs[3] * rhs.coeffs[1];
        coeffs[3] = tmp.coeffs[0] * rhs.coeffs[3] + tmp.coeffs[3] * rhs.coeffs[0];

        return *this;
    }

    E2MV<T>& operator^=(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;

        coeffs[3] = tmp.coeffs[1] * rhs.coeffs[2] - tmp.coeffs[2] * rhs.coeffs[1];

        return *this;
    }

    E2MV<T>& operator*=(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;

        coeffs[0] = tmp.coeffs[0] * rhs.coeffs[0] + tmp.coeffs[1] * rhs.coeffs[1] + tmp.coeffs[2] * rhs.coeffs[2] - tmp.coeffs[3] * rhs.coeffs[3];
        coeffs[1] = tmp.coeffs[0] * rhs.coeffs[1] + tmp.coeffs[1] * rhs.coeffs[0] - tmp.coeffs[2] * rhs.coeffs[3] + tmp.coeffs[3] * rhs.coeffs[2];
        coeffs[2] = tmp.coeffs[0] * rhs.coeffs[2] + tmp.coeffs[1] * rhs.coeffs[3] + tmp.coeffs[2] * rhs.coeffs[0] - tmp.coeffs[3] * rhs.coeffs[1];
        coeffs[3] = tmp.coeffs[0] * rhs.coeffs[3] + tmp.coeffs[3] * rhs.coeffs[0] + tmp.coeffs[1] * rhs.coeffs[2] - tmp.coeffs[2] * rhs.coeffs[1];

        return *this;
    }

};