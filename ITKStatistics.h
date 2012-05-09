#ifndef ITKStatistics_H
#define ITKStatistics_H

namespace ITKStatistics
{

template<typename TVector>
typename TypeTraits<TVector>::LargerComponentType RunningAverage(const TVector& v);

template<typename TVector>
typename TypeTraits<TVector>::LargerComponentType Average(const TVector& v);

template<typename TVector>
typename TypeTraits<TVector>::LargerComponentType Variance(const TVector& v);

}

#include "ITKStatistics.hpp"

#endif
