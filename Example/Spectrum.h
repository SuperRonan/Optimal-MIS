#pragma once

namespace MISExample
{

    // Example spectrum class used by the optimal MIS estimator
    // For the different estimator class to compile, 
    // the Spectrum template argument should provide at least the functionalities of this implementation
    template<class Float, int N>
    class Spectrum
    {
    protected:

#define me (*this)

        Float m_data[N];

    public:

        // This one is necessary for the estimator to know how much memory to allocate
        // So the static constexpr are also very important
        static constexpr int size()
        {
            return N;
        }

        Spectrum(Float f=0)
        {
            std::fill_n(m_data, N, f); 
        }

        Spectrum(Spectrum const& other)
        {
            std::copy(other.m_data, other.m_data+N, m_data);
        }

        Spectrum& operator=(Spectrum const& other)
        {
            std::copy(other.m_data, other.m_data + N, m_data);
            return me;
        }

        bool operator==(Spectrum const& other)const
        {
            for (int i = 0; i < N; ++i)
                if (me[i] != other[i])   return false;
            return true;
        }

        bool operator!=(Spectrum const& other)const
        {
            return !(me == other);
        }

        bool isZero()const
        {
            return me == Float(0);
        }


        Float& operator[](int i)
        {
            return m_data[i];
        }

        Float const& operator[](int i)const
        {
            return m_data[i];
        }

        // basic mathematical operators

        Spectrum& operator+=(Spectrum const& other)
        {
            for (int i = 0; i < N; ++i)
                me[i] += other[i];
            return me;
        }

        Spectrum& operator-=(Spectrum const& other)
        {
            for (int i = 0; i < N; ++i)
                me[i] -= other[i];
            return me;
        }

        Spectrum& operator*=(Spectrum const& other)
        {
            for (int i = 0; i < N; ++i)
                me[i] *= other[i];
            return me;
        }

        Spectrum& operator/=(Spectrum const& other)
        {
            for (int i = 0; i < N; ++i)
                me[i] /= other[i];
            return me;
        }

        Spectrum& operator*=(Float f)
        {
            for (int i = 0; i < N; ++i)
                me[i] *= f;
            return me;
        }

        Spectrum& operator/=(Float f)
        {
            for (int i = 0; i < N; ++i)
                me[i] /= f;
            return me;
        }

        Spectrum operator+(Spectrum const& other)const
        {
            Spectrum res = me;
            res += other;
            return res;
        }

        Spectrum operator-(Spectrum const& other)const
        {
            Spectrum res = me;
            res -= other;
            return res;
        }

        Spectrum operator*(Spectrum const& other)const
        {
            Spectrum res = me;
            res *= other;
            return res;
        }

        Spectrum operator/(Spectrum const& other)const
        {
            Spectrum res = me;
            res /= other;
            return res;
        }

        Spectrum operator*(Float f)const
        {
            Spectrum res = me;
            res *= f;
            return res;
        }

        Spectrum operator/(Float f)const
        {
            Spectrum res = me;
            res /= f;
            return res;
        }
        
#undef me
    };
}


template <class Stream, class Float, int N>
Stream& operator<<(Stream& stream, MISExample::Spectrum<Float, N> const& spec)
{
    stream << '[';
    for (int i = 0; i < spec.size(); ++i)
    {
        stream << std::to_string(spec[i]);
        if (i < spec.size() - 1)
            stream << ", ";
    }
    stream << "]";
    return stream;
}