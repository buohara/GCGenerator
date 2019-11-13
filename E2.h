static const uint32_t numCoeffs = 4;

template<typename T> struct E2MV
{
    T coeffs[numCoeffs];

    E2MV<T>() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }

    E2MV<T>(T S, T E0, T E1, T E0E1)
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

    E2MV<T>& operator*=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] *= rhs;

        return *this;
    }

    E2MV<T>& operator/=(const E2MV<T>& rhs)
    {
        E2MV<T>tmp = ~rhs;
        T mag = rhs.Mag();
        (*this) *= tmp;
        (*this) /= (mag * mag);
        return *this;
    }

    E2MV<T>& operator/=(const T rhs)
    {
        for (uint32_t i = 0; i < numCoeffs; i++)
            coeffs[i] /= rhs;

        return *this;
    }

    E2MV<T> operator+(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E2MV<T> operator-(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;
        tmp -= rhs;
        return tmp;
    }

    E2MV<T> operator|(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;
        tmp |= rhs;
        return tmp;
    }

    E2MV<T> operator^(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;
        tmp ^= rhs;
        return tmp;
    }

    E2MV<T> operator*(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;
        tmp *= rhs;
        return tmp;
    }

    E2MV<T> operator*(const T rhs)
    {
        E2MV<T> tmp = *this;
        tmp += rhs;
        return tmp;
    }

    E2MV<T> operator/(const E2MV<T>& rhs)
    {
        E2MV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E2MV<T> operator/(const T rhs)
    {
        E2MV<T> tmp = *this;
        tmp /= rhs;
        return tmp;
    }

    E2MV<T> operator~()
    {
        E2MV<T> tmp = *this;
        tmp.coeffs[3] = -tmp.coeffs[3];

        return tmp;
    }

    T S(){ return coeffs[0]; }
    void S(T val){ coeffs[0] = val; }

    T E0(){ return coeffs[1]; }
    void E0(T val){ coeffs[1] = val; }

    T E1(){ return coeffs[2]; }
    void E1(T val){ coeffs[2] = val; }

    T E0E1(){ return coeffs[3]; }
    void E0E1(T val){ coeffs[3] = val; }

    T Mag()
    {
        T mag = 0.0;
        for (uint32_t i = 0; i < nCoeffs; i++) mag += coeffs[i] * coeffs[i];
        return sqrt(mag);
    }

};