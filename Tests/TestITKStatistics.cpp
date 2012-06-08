// NOTE! ITKContainerInterface must be included BEFORE Statistics or it will not know about the
// definitions of the functions required and produce compiler errors.
#include "ITKContainerInterface.h"
#include "Helpers/Statistics.h"

// ITK
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
