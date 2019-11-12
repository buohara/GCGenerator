static const uint32_t numCoeffs = 4;

template<typename T> struct E2MV
{
    T coeffs[4];

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

};