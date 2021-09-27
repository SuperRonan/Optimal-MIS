#pragma once

#include <type_traits>
#include <cassert>

namespace MIS
{
	/// <summary>
	/// Wrapper class for the Estimators. Allows to use undefined function (operator[], size()) on double/float.
	/// </summary>
	/// <typeparam name="Spectrum"></typeparam>
	template <class Spectrum>
	class SpectrumWrapper
	{
	protected:

		static constexpr bool isUnique() // Intended for double of float
		{
			return std::is_arithmetic<Spectrum>::value;
		}

		Spectrum & m_data;

	public:

		SpectrumWrapper(Spectrum & data):
			m_data(data)
		{}

		constexpr static int size()
		{
			if constexpr (isUnique())
			{
				return 1;
			}
			else
				return Spectrum::size();
		}

		operator Spectrum&()const
		{
			return m_data;
		}

		auto& operator[](int i)
		{
			if constexpr (isUnique())
			{
				assert(i == 0);
				return m_data;
			}
			else
				return m_data[i];
		}

		auto operator[](int i)const
		{
			if constexpr (isUnique())
			{
				assert(i == 0);
				return m_data;
			}
			else
				return m_data[i];
		}
		
		bool isZero()const
		{
			return m_data == 0;
		}
	};
}