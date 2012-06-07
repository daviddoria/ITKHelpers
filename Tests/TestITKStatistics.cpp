#include "ITKStatistics.h"

#include "itkVariableLengthVector.h"
#include "itkCovariantVector.h"

int main()
{
  {
  itk::VariableLengthVector<float> v(5);
  std::cout << Statistics::Average(v);
  std::cout << Statistics::RunningAverage(v);
  std::cout << Statistics::Variance(v);
  }
  
  {
  itk::CovariantVector<float, 3> v;
  std::cout << Statistics::Average(v);
  std::cout << Statistics::RunningAverage(v);
  std::cout << Statistics::Variance(v);
  }
  
  return 0;
}
