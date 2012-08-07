/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef ITKHelpers_H
#define ITKHelpers_H

#include "ITKHelpersTypes.h"

// Submodules
class Mask;
#include <Helpers/TypeTraits.h>

// STL
#include <string>

// ITK
#include "itkImage.h"
#include "itkImageIOBase.h" // For GetPixelTypeFromFile
#include "itkIndex.h"
#include "itkRGBPixel.h"
#include "itkSize.h"
#include "itkVectorImage.h"

namespace ITKHelpers
{
  using namespace ITKHelpersTypes;
////////////////////////////////////////////////////////////////////////
///////// Function templates (defined in ITKHelpers.hpp) /////////
////////////////////////////////////////////////////////////////////////

/** Convert an RGB image to the CIELAB colorspace.
 * 'rgbImage' cannot be const because the adaptor doesn't allow it. */
template<typename TOutputImage>
void RGBImageToHSVImage(RGBImageType* const rgbImage, TOutputImage* const hsvImage);

/** Convert an image to the HSV colorspace. */
template <typename TInputImage, typename TOutputImage>
void ITKImageToHSVImage(const TInputImage* const image, TOutputImage* const hsvImage);

/** Convert an RGB image to the CIELAB colorspace.
 * 'rgbImage' cannot be const because the adaptor doesn't allow it. */
template<typename TOutputImage>
void RGBImageToCIELabImage(RGBImageType* const rgbImage, TOutputImage* const cielabImage);

/** Determine if any of the 8 neighbors pixels has the specified value. */
template<typename TImage>
bool HasNeighborWithValue(const itk::Index<2>& pixel, const TImage* const image,
                          const typename TImage::PixelType& value);

/** Blur all channels of an image. */
template<typename TImage>
void BlurAllChannels(const TImage* const image, TImage* const output,
                     const float sigma);

/** Blur all channels of an image, preserving edges. */
template<typename TInputImage, typename TPixelType>
void BilateralFilterAllChannels(const TInputImage* const image, itk::VectorImage<TPixelType,2>* const output,
                                const float domainSigma, const float rangeSigma);

/** Set the values of the pixels on the boundary of the 'region' to 'value'. */
template<typename TImage>
void OutlineRegion(TImage* image, const itk::ImageRegion<2>& region, const typename TImage::PixelType& value);

/** Deep copy a scalar image. */
template<typename TInputPixel, typename TOutputPixel>
void DeepCopy(const itk::Image<TInputPixel, 2>* const input, itk::Image<TOutputPixel, 2>* const output);

/** Deep copy a vector image. */
template<typename TInputPixel, typename TOutputPixel>
void DeepCopy(const itk::VectorImage<TInputPixel, 2>* const input,
              itk::VectorImage<TOutputPixel, 2>* const output);

/** Deep copy a 'region' of 'input' to the same 'region' in 'output'. */
template<typename TInputImage, typename TOutputImage>
void DeepCopyInRegion(const TInputImage* const input, const itk::ImageRegion<2>& region,
                      TOutputImage* const output);

/** Copy the region around 'sourcePosition' with radius 'radius' to the region around
  * 'targetPosition' with radius 'radius'. */
template <class TImage>
void CopyRegion(const TImage* const sourceImage, TImage* const targetImage,
                const itk::Index<2>& sourcePosition,
                const itk::Index<2>& targetPosition, const unsigned int radius);

/** Copy 'sourceRegion' in 'sourceImage' to 'targetRegion' in 'targetImage'. */
template <class TImage>
void CopyRegion(const TImage* const sourceImage, TImage* const targetImage, const itk::ImageRegion<2>& sourceRegion,
               const itk::ImageRegion<2>& targetRegion);

/** The same as the other function by the same name, but with the argument order switched. */
template <class TImage>
void CopyRegion(const TImage* const sourceImage, const itk::ImageRegion<2>& sourceRegion,
                TImage* const targetImage, const itk::ImageRegion<2>& targetRegion);

/** Copy 'sourceRegion' to 'targetRegion' in 'image'. */
template <class TImage>
void CopySelfRegion(TImage* const image, const itk::ImageRegion<2>& sourceRegion,
                    const itk::ImageRegion<2>& targetRegion);

/** Replace pixels in 'image' with value 'queryValue' with 'replacementValue'. */
template<typename TImage>
void ReplaceValue(TImage* image, const typename TImage::PixelType& queryValue,
                  const typename TImage::PixelType& replacementValue);

/** Find the min value of 'image' in 'region' of a scalar image. */
template <class TImage>
typename TImage::PixelType MinValue(const TImage* const image, const itk::ImageRegion<2>& region);

/** Find the max value of 'image' in 'region' of a scalar image. */
template <class TImage>
typename TImage::PixelType MaxValue(const TImage* const image, const itk::ImageRegion<2>& region);

/** Find the min value of 'image' of a scalar image. */
template <class TImage>
typename TImage::PixelType MinValue(const TImage* const image);

/** Find the max value of 'image' of a scalar image. */
template <class TImage>
typename TImage::PixelType MaxValue(const TImage* const image);

/** Get the maximum values of every channel of a vector image. */
template <class TPixel>
std::vector<TPixel> MaxValues(const itk::VectorImage<TPixel, 2>* const image);

/** Compute the max value of each channel of 'image' in 'region'. */
template <class TPixel>
std::vector<TPixel> MaxValues(const itk::VectorImage<TPixel, 2>* const image, const itk::ImageRegion<2>& region);

/** Compute the min value of each channel of 'image' in 'region'. */
template <class TPixel>
std::vector<TPixel> MinValues(const itk::VectorImage<TPixel, 2>* const image, const itk::ImageRegion<2>& region);

/** Compute the min value of each channel of the image. */
template <class TPixel>
std::vector<TPixel> MinValues(const itk::VectorImage<TPixel, 2>* const image);

/** Find the location of the pixel with the maximum value. */
template <class TImage>
itk::Index<2> MaxValueLocation(const TImage* const image);

/** Find the location of the pixel with the minimum value. */
template <class TImage>
itk::Index<2> MinValueLocation(const TImage* const image);

/** Convert a color image to a grayscale image. */
template <typename TInputImage, typename TOutputImage>
void ColorToGrayscale(const TInputImage* const colorImage, TOutputImage* const grayscaleImage);

/** Set the pixels inside 'region' in 'image' to 'blankValue', and those on the border of the 'region'
 *  to 'outlineValue'. */
template<typename TImage>
void BlankAndOutlineRegion(TImage* const image, const itk::ImageRegion<2>& region,
                           const typename TImage::PixelType& blankValue,
                           const typename TImage::PixelType& outlineValue);

template<typename TImage>
void SetRegionToConstant(TImage* const image, const itk::ImageRegion<2>& region,
                         const typename TImage::PixelType& constant);

template<typename TImage>
void SetImageToConstant(TImage* const image, const typename TImage::PixelType& constant);

template<typename TImage>
unsigned int CountNonZeroPixels(const TImage* const image);

template<typename TImage>
unsigned int CountPixelsWithValue(const TImage* const image, const typename TImage::PixelType& value);

/** Get the non-zero pixels in an image. */
template<typename TImage>
std::vector<itk::Index<2> > GetNonZeroPixels(const TImage* const image);

/** Get the non-zero pixels in a region. AKA BinaryImageToPixelList */
template<typename TImage>
std::vector<itk::Index<2> > GetNonZeroPixels(const TImage* const image, const itk::ImageRegion<2>& region);

/** Set the region of 'image' and initialize it to zero. */
template<typename TImage>
void InitializeImage(TImage* const image, const itk::ImageRegion<2>& region);

/** Initialize a vector image to zero. */
template<typename TPixel>
void InitializeImage(const itk::VectorImage<TPixel, 2>* const input, const unsigned int numberOfComponents,
                     const itk::ImageRegion<2>& region);

/** Dilate 'image' in-place. */
template<typename TImage>
void DilateImage(const TImage* const image, const unsigned int radius);

/** Dilate 'image' and save the result in 'output'. */
template<typename TImage>
void DilateImage(const TImage* const image, const unsigned int radius, TImage* const dilatedImage);

/** Subtract regions. */
template<typename TImage>
void SubtractRegions(const TImage* const image1, const itk::ImageRegion<2>& region1,
                     const TImage* const image2, const itk::ImageRegion<2>& region2, TImage* const output);

/** AND regions. */
template<typename TImage>
void ANDRegions(const TImage* const image1, const itk::ImageRegion<2>& region1,
                const TImage* const image2, const itk::ImageRegion<2>& region2, itk::Image<bool, 2>* const output);

/** XOR images. */
template<typename TInputImage, typename TOutputImage>
void XORImages(const TInputImage* const image1, const TInputImage* const image2, TOutputImage* const output,
               const typename TOutputImage::PixelType& trueValue = 1);

/** XOR regions. */
template<typename TInputImage, typename TOutputImage>
void XORRegions(const TInputImage* const image1, const itk::ImageRegion<2>& region1,
                const TInputImage* const image2, const itk::ImageRegion<2>& region2, TOutputImage* const output,
                const typename TOutputImage::PixelType& trueValue = 1);

/** Change the value of all pixels with value = 'oldValue' to 'newValue. */
template<typename TImage>
void ChangeValue(const TImage* const image, const typename TImage::PixelType& oldValue,
                 const typename TImage::PixelType& newValue);

/** Extract a channel of an image. The output image should be a scalar image,
  * but does not have to have the same pixel type as the input image. */
template<typename TInputImage, typename TOutputImage>
void ExtractChannel(const TInputImage* const image, const unsigned int channel,
                    TOutputImage* const output);

/** A specialization to allow this function to pass through a scalar image. */
template<typename TInputPixel, typename TOutputPixel>
void ExtractChannel(const itk::Image<TInputPixel, 2>* const image, const unsigned int channel,
                    itk::Image<TOutputPixel, 2>* const output);

/** An even further specialization to allow this function to process Image<CovariantVector<> >. */
template<typename TInputPixelComponent, unsigned int PixelDimension, typename TOutputPixel>
void ExtractChannel(const itk::Image<itk::CovariantVector<TInputPixelComponent, PixelDimension>, 2>* const image, const unsigned int channel,
                    itk::Image<TOutputPixel, 2>* const output);

/** An even further  specialization to allow this function to process Image<Vector<> >. */
template<typename TInputPixelComponent, unsigned int PixelDimension, typename TOutputPixel>
void ExtractChannel(const itk::Image<itk::Vector<TInputPixelComponent, PixelDimension>, 2>* const image, const unsigned int channel,
                    itk::Image<TOutputPixel, 2>* const output);

/** Extract a channels of an image. */
template<typename TInputImage, typename TOutputImage>
void ExtractChannels(const TInputImage* const image, const std::vector<unsigned int> channels,
                    TOutputImage* const output);

/** Extract a region of an image. */
template<typename TImage>
void ExtractRegion(const TImage* const image, const itk::ImageRegion<2>& region,
                   TImage* const output);

/** Scale a channel of an image between 0 and 'channelMax'. */
template<typename TPixel>
void ScaleChannel(const itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                  const TPixel channelMax, typename itk::VectorImage<TPixel, 2>* const output);

/** Force an image to be 3 channels. If it is already 3 channels, copy it through. If it is less than 3 channels. */
template<typename TPixel>
void ConvertTo3Channel(const typename itk::VectorImage<TPixel, 2>* const image,
                      typename itk::VectorImage<TPixel, 2>* const output);

/** Replace a channel of an image. */
template<typename TPixel>
void ReplaceChannel(const typename itk::VectorImage<TPixel, 2>* const image, const unsigned int channel,
                    const typename itk::Image<TPixel, 2>* const replacement,
                    typename itk::VectorImage<TPixel, 2>* const output);

/** A specialization to allow this function to pass through scalar images. */
template<typename TPixel>
void ReplaceChannel(const itk::Image<TPixel, 2>* const image, const unsigned int channel,
                    const itk::Image<TPixel, 2>* const replacement,
                    itk::Image<TPixel, 2>* const output);

/** Read an image from a file. */
template<typename TImage>
void ReadImage(const std::string&, TImage* const image);

/** Normalize every pixel/vector in a vector image. */
template<typename TImage>
void NormalizeVectorImage(TImage* const image);

/** Create an itk::Index<2> from any object with a operator[]. */
template<typename T>
itk::Index<2> CreateIndex(const T& v);

/** Get the average value of the neighbors of a pixel. */
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType AverageNeighborValue(const TImage* const image,
                                                                                 const itk::Index<2>& pixel);

/** Get the locations of the 8-neighors of 'pixel' with value 'value'. */
template<typename TImage>
std::vector<itk::Index<2> > Get8NeighborsWithValue(const itk::Index<2>& pixel, const TImage* const image,
                                                   const typename TImage::PixelType& value);

/** Get the locations of the 8-neighors of 'pixel' with value 'value', only if they are inside 'region'. */
template<typename TImage>
std::vector<itk::Index<2> > Get8NeighborsInRegionWithValue(const itk::Index<2>& pixel, const TImage* const image,
                                                           const itk::ImageRegion<2>& region,
                                                           const typename TImage::PixelType& value);

/** Get the locations of the 8-neighors of 'pixel' that do NOT equal 'value', only if they are inside 'region'. */
template<typename TImage>
std::vector<itk::Index<2> > Get8NeighborsInRegionNotEqualValue(const itk::Index<2>& pixel, const TImage* const image,
                                                               const itk::ImageRegion<2>& region,
                                                               const typename TImage::PixelType& value);

/** Get the locations of pixels with value 'value' in 'region'. */
template<typename TImage>
std::vector<itk::Index<2> > GetPixelsWithValue(const TImage* const image, const itk::ImageRegion<2>& region,
                                               const typename TImage::PixelType& value);

/** Get the locations of pixels with value 'value'. */
template<typename TImage>
std::vector<itk::Index<2> > GetPixelsWithValue(const TImage* const image,
                                               const typename TImage::PixelType& value);

/** Fetch the values in the image at the specified indices. AKA GetPixelsAtIndices */
template<typename TImage>
std::vector<typename TImage::PixelType> GetPixelValues(const TImage* const image,
                                                       const std::vector<itk::Index<2> >& indices);

/** Linearize the values in the image region. */
template<typename TImage>
std::vector<typename TImage::PixelType> GetPixelValuesInRegion(const TImage* const image,
                                                               const itk::ImageRegion<2>& region);

/** Compute the average of the values appearing at the specified indices. */
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType AverageOfPixelsAtIndices(
  const TImage* const image, const std::vector<itk::Index<2> >& indices);

/** Compute the variance of the values appearing at the specified indices. The variance of the ith component is the
 * ith component of the output pixel*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType
VarianceOfPixelsAtIndices(const TImage* const image, const std::vector<itk::Index<2> >& indices);

/** Compute the average of all pixels in an image.*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType
AverageOfImage(const TImage* const image);

/** Compute the average of all pixels in a region.*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType
AverageInRegion(const TImage* const image, const itk::ImageRegion<2>& region);

/** Compute the variance of all pixels in an image.*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType VarianceOfImage(const TImage* const image);

/** Compute the variance of all pixels in a region.*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType
VarianceInRegion(const TImage* const image, const itk::ImageRegion<2>& region);

/** Compute the average difference of corresponding pixels from two regions of an image.*/
template<typename TImage, typename TDifferenceFunctor>
float AverageDifferenceInRegion(const TImage* const image, const itk::ImageRegion<2>& region1,
                                const itk::ImageRegion<2>& region2, TDifferenceFunctor differenceFunctor);

/** Compute the average difference of corresponding pixels from regions in two images.*/
template<typename TImage, typename TDifferenceFunctor>
float AverageDifferenceInRegion(const TImage* const image1, const itk::ImageRegion<2>& region1,
                                const TImage* const image2, const itk::ImageRegion<2>& region2,
                                TDifferenceFunctor differenceFunctor);

/** Sum the componets of an object.*/
template <typename T>
typename TypeTraits<T>::LargerComponentType SumOfComponents(const T& v);

/** Sum the magnitudes of the componets of an object.*/
template <typename T>
typename TypeTraits<T>::LargerComponentType SumOfComponentMagnitudes(const T& v);

/** Set any object to all zeros. Usually scalars or vectors */
template<typename T>
void SetObjectToZero(T& object);

/** Print the image on the screen */
template<typename TImage>
void PrintImage(const TImage* const image);

/** Print a region of an image on the screen */
template<typename TImage>
void PrintRegion(const TImage* const image, const itk::ImageRegion<2>& region);

template<typename TVectorImage, typename TScalarImage>
void SetChannel(TVectorImage* const vectorImage, const unsigned int channel, const TScalarImage* const image);

template<typename TImage>
void ScaleAllChannelsTo255(TImage* const image);

template<typename TInputImage,typename TOutputImage>
void CastImage(const TInputImage* const inputImage, TOutputImage* const outputImage);

template<typename TImage>
void ScaleTo255(TImage* const image);

/** Set all of the 'pixels' to 'value'. AKA PixelList */
template <typename TImage, typename TPixel>
void SetPixels(TImage* const image, const std::vector<itk::Index<2> >& pixels, const TPixel& value);

/** Set all pixels in 'region' to 'value' */
template <typename TImage>
void SetPixelsInRegionToValue(TImage* const image, const itk::ImageRegion<2>& region,
                              const typename TImage::PixelType& value);

/** Extract a region and then normalize it. */
template<typename TImage>
void ExtractAndNormalizeRegion(const TImage* const image, const itk::ImageRegion<2>& region,
                               TImage* const outputImage);

/** For each component of each pixel, subtract the mean of the channel and divide by
  * the standard deviation of the corresponding channel. */
template<typename TImage>
void NormalizeImageChannels(const TImage* const image, TImage* const outputImage);

/** Compute the mean of the image. */
template<typename TImage>
float MeanValue(const TImage* const image);

/** Compute the standard deviation of the image. */
template<typename TImage>
float StandardDeviation(const TImage* const image);

/** Compute the variance of the image. */
template<typename TImage>
float Variance(const TImage* const image);

/** Compute the minimum value of each channel of the image and return them in a vector whose length is
 *  the same as the number of components of the image. */
template<typename TImage>
std::vector<typename TImage::InternalPixelType> ComputeMinOfAllChannels(const TImage* const image);

/** Compute the maximum value of each channel of the image and return them in a vector whose length is
 *  the same as the number of components of the image. */
template<typename TImage>
std::vector<typename TImage::InternalPixelType> ComputeMaxOfAllChannels(const TImage* const image);

/** Create a string representation of a vector. */
template <typename TVector>
std::string VectorToString(const TVector& vec);

/** Create a string representation of an itk::CovariantVector. */
template <typename TComponent, unsigned int Dimension>
std::string VectorToString(const itk::CovariantVector<TComponent, Dimension>& vec);

/** Downsample 'image' by 'factor'. */
template <typename TImage>
void Downsample(const TImage* const image, const float factor, TImage* const output);

/** Upsample 'image' by 'factor'. */
template <typename TImage>
void Upsample(const TImage* const image, const float factor, TImage* const output);

/** Strip a row and/or column to make both dimensions even. */
template <typename TImage>
void CreateEvenSizeImage(const TImage* const image, TImage* const output);

/** Scale 'image' to match the 'destinationSize'. */
template <typename TImage>
void ScaleImage(const TImage* const image, const itk::Size<2>& destinationSize, TImage* const output);

template <typename TInputImage, typename TOutputImage>
void MagnitudeImage_Generic(const TInputImage* const image, TOutputImage* const output);

/** Compute the magnitude image. */
template <typename TInputImage, typename TOutputImage>
void MagnitudeImage(const TInputImage* const image, TOutputImage* const output);

/** If MagnitudeImage is called on a scalar image, compute the absolute value rather than the magnitude. */
template <typename TInputPixel, typename TOutputPixel>
void MagnitudeImage(const itk::Image<TInputPixel, 2>* const image, itk::Image<TOutputPixel,2>* const output);

/** Compute the magnitude image. */
template <typename TInputPixel, unsigned int TVectorDim, typename TOutputPixel>
void MagnitudeImage(const itk::Image<itk::CovariantVector<TInputPixel, TVectorDim>, 2>* const image,
                    itk::Image<TOutputPixel,2>* const output);

/**  Write the first 3 channels of 'image' as an unsigned char (RGB) image. */
template <typename TPixel>
void WriteVectorImageAsRGB(const itk::VectorImage<TPixel,2>* const image, const std::string& fileName);

/** Write a 'region' of an 'image' to 'filename'.*/
template <typename TPixel>
void WriteVectorImageRegionAsRGB(const itk::VectorImage<TPixel,2>* const image,
                                 const itk::ImageRegion<2>& region, const std::string& filename);

/** Compute a histogram of gradients. */
template <typename TImage>
std::vector<float> HistogramOfGradients(const TImage* const image,
                                        const itk::ImageRegion<2>& region, const unsigned int numberOfBins);

/** Compute the gradient of an image. */
template <typename TImage>
void ComputeGradients(const TImage* const image,
                      itk::Image<itk::CovariantVector<float, 2>, 2>* output);

/** Compute a histogram of gradients when the gradients are already known. */
template <typename TGradientImage>
std::vector<float> HistogramOfGradientsPrecomputed(const TGradientImage* const image,
                                                   const itk::ImageRegion<2>& region,
                                                   const unsigned int numberOfBins);

/** Convert an RGB image to a 3-channel VectorImage. */
template <typename TOutputPixel>
void RGBImageToVectorImage(const itk::Image<itk::RGBPixel<unsigned char>, 2>* const image,
                           itk::VectorImage<TOutputPixel, 2>* const outputImage);

/** Convert the first 3 channels of an ITK image to the CIELAB colorspace. */
template <typename TInputImage, typename TOutputImage>
void ITKImageToCIELabImage(const TInputImage* const rgbImage, TOutputImage* const cielabImage);

/** Compute the central derivative of 'image' in 'direction' and store it in 'output'. */
template <class TImage>
void CentralDifferenceDerivative(const TImage* const image, const unsigned int direction, TImage* const output);

/** Find the first non-zero pixel in raster scan order. */
template <class TImage>
itk::Index<2> FindFirstNonZeroPixel(const TImage* const image);

/** Get open contour ordering. */
template <class TImage>
std::vector<itk::Index<2> > GetOpenContourOrdering(const TImage* const image, const itk::Index<2>& start);

/** Get closed contour ordering. */
template <class TImage>
std::vector<itk::Index<2> > GetClosedContourOrdering(const TImage* const image, const itk::Index<2>& start);

/** Construct a breadth first (BFS) ordering on the non-zero pixels attached to a 'start' pixel. */
template <class TImage>
std::vector<itk::Index<2> > BreadthFirstOrderingNonZeroPixels(const TImage* const image, const itk::Index<2>& start);

/** Check if the non-zero pixels attached to a 'start' pixel form a closed loop. */
template <class TImage>
bool IsClosedLoop(const TImage* const image, const itk::Index<2>& start);

/**  Write the first 3 channels of an image to a file as unsigned chars. */
template<typename TImage>
void WriteRGBImage(const TImage* const input, const std::string& filename);

/** Write an image to a file named 'prefix'_'iteration'.extension*/
template <typename TImage>
void WriteSequentialImage(const TImage* const image, const std::string& filePrefix, const unsigned int iteration,
                          const unsigned int iterationLength, const std::string& extension);

/** Write 'image' to 'fileName' if 'condition' is true. */
template <typename TImage>
void WriteImageConditional(const TImage* const image, const std::string& fileName, const bool condition);

/** Scale a scalar image to the range (0-255) and then write it to a file as an unsigned char image. */
template <class TImage>
void WriteScaledScalarImage(const TImage* const image, const std::string& filename);

/** Write 'image' to 'fileName'.*/
template<typename TImage>
void WriteImage(const TImage* const image, const std::string& fileName);

/** Write a 'region' of an 'image' to 'filename'.*/
template<typename TImage>
void WriteRegion(const TImage* const image, const itk::ImageRegion<2>& region, const std::string& filename);

/** Write a vector 'v' to the screen.*/
template <typename T>
void OutputVector(const std::vector<T>& v);

/** Combine two VectorImages into a VectorImage where the number of channels is
  * the sum of the two input images. */
template <typename TPixel>
void StackImages(const typename itk::VectorImage<TPixel, 2>* const image1,
                 const typename itk::VectorImage<TPixel, 2>* const image2,
                 typename itk::VectorImage<TPixel, 2>* const output);

/** Convert the first 3 channels of a float vector image to an unsigned char/color/rgb image. */
template <typename TImage>
void VectorImageToRGBImage(const TImage* const image, RGBImageType* const rgbImage);

/** Draw a rectangle in an image. */
template <typename TImage>
void DrawRectangle(TImage* const image, const typename TImage::PixelType& value,
                   const itk::Index<2>& corner0, const itk::Index<2>& corner1);

/** Interpolate points between p0 and p1 by weighting the endpoint values by their
  * distance to the current pixel.*/
template<typename TImage>
void InterpolateLineBetweenPoints(TImage* const image, const itk::Index<2>& p0, const itk::Index<2>& p1);

/** Generate a random image. */
template<typename TPixel>
void RandomImage(itk::Image<TPixel, 2>* const image);

template<typename TPixel>
void RandomImage(itk::VectorImage<TPixel, 2>* const image);

/** Find which element in 'vec' is closest to 'value'. */
template<typename TValue>
unsigned int ClosestValueIndex(const std::vector<TValue>& vec,
                               const TValue& value);

/** Find which element in 'vec' is closest to 'value'. Specialized for
  * itk::CovariantVector (because it has a special GetSquaredNorm() function). */
template<typename TComponent, unsigned int NumberOfComponents>
unsigned int ClosestValueIndex(
   const std::vector<itk::CovariantVector<TComponent, NumberOfComponents> >& vec,
   const itk::CovariantVector<TComponent, NumberOfComponents>& value);

/** Compute the bounding box of all pixels with 'value'. */
template<typename TImage>
itk::ImageRegion<2> ComputeBoundingBox(const TImage* const image,
                                       const typename TImage::PixelType& value);


//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Non-template function declarations (defined in Helpers.cpp) ///////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

/** Mark each pixel at the specified 'indices' as a non-zero pixel in 'image' */
void IndicesToBinaryImage(const std::vector<itk::Index<2> >& indices,
                          UnsignedCharScalarImageType* const image);

/** Get the number of components per pixel in an image file. */
unsigned int GetNumberOfComponentsPerPixelInFile(const std::string& filename);

/** Get a short string of an itk::Index */
std::string GetIndexString(const itk::Index<2>& index);

/** Get a short string of an itk::Size */
std::string GetSizeString(const itk::Size<2>& size);

/** Get a region from the corner (0,0) with size 'size'. */
itk::ImageRegion<2> CornerRegion(const itk::Size<2>& size);

/** If an image is one of our few predefined types, print it's type the the screen. */
void OutputImageType(const itk::ImageBase<2>* const input);

/** Average each component of a list of vectors then construct and return a new vector
  * composed of these averaged components. */
FloatVector2Type AverageVectors(const std::vector<FloatVector2Type>& vectors);

/** This function creates an Offset<2> by setting the 'dimension' index of the Offset<2>
 * to the value of the Offset<1>, and filling the other position with a 0. */
itk::Offset<2> OffsetFrom1DOffset(const itk::Offset<1>& offset1D, const unsigned int dimension);

/** Compute the angle between two vectors. */
float AngleBetween(const FloatVector2Type& v1, const FloatVector2Type& v2);

/** Get the 'quadrant'th quadrant (quarter) region of the 'region'. */
itk::ImageRegion<2> GetQuadrant(const itk::ImageRegion<2>& region, const unsigned int quadrant);

/** Get the center pixel of a region. The region is assumed to have odd dimensions. */
itk::Index<2> GetRegionCenter(const itk::ImageRegion<2>& region);

/** Get the size of a region by it's radius. E.g. radius 5 -> 11x11 patch. */
itk::Size<2> SizeFromRadius(const unsigned int radius);

/** "Follow" a vector from one pixel to find the next pixel it would "hit". */
itk::Index<2> GetNextPixelAlongVector(const itk::Index<2>& pixel, const FloatVector2Type& vector);

/** Get the direction in integer pixel coordinates of the 'vector' */
itk::Offset<2> GetOffsetAlongVector(const FloatVector2Type& vector);

/** Make an ImageRegion centered on 'pixel' with radius 'radius'. */
itk::ImageRegion<2> GetRegionInRadiusAroundPixel(const itk::Index<2>& pixel, const unsigned int radius);

/** Get the offsets of the 8 neighborhood of a pixel. */
std::vector<itk::Offset<2> > Get8NeighborOffsets();

/** Get the indices of the neighbors of a 'pixel'. */
std::vector<itk::Index<2> > Get8Neighbors(const itk::Index<2>& pixel);

/** Get the indices of the neighbors of a 'pixel' that are inside of a 'region'. */
std::vector<itk::Index<2> > Get8NeighborsInRegion(const itk::ImageRegion<2>& region, const itk::Index<2>& pixel);

/** Get the neighboring of a 'pixel' (of size 'queryRegionSize') that are inside of a 'searchRegion'. */
std::vector<itk::ImageRegion<2> > Get8NeighborRegionsInRegion(const itk::ImageRegion<2>& searchRegion, const itk::Index<2>& pixel,
                                                        const itk::Size<2>& queryRegionSize);

/** The return value MUST be a smart pointer. */
itk::ImageBase<2>::Pointer CreateImageWithSameType(const itk::ImageBase<2>* input);

/** Convert a list of offsets from 'index' to a list of absolute indices. */
std::vector<itk::Index<2> > OffsetsToIndices(const std::vector<itk::Offset<2> >& offsets,
                                             const itk::Index<2>& index);

/** Interpret a list of offsets as a list of absolute indices. */
std::vector<itk::Index<2> > OffsetsToIndices(const std::vector<itk::Offset<2> >& offsets);

/** Interpret a list of indices as a list of offsets. */
std::vector<itk::Offset<2> > IndicesToOffsets(const std::vector<itk::Index<2> >& indices,
                                              const itk::Index<2>& index);

/** Get the locations of the pixels on the boundary of a region. */
std::vector<itk::Index<2> > GetBoundaryPixels(const itk::ImageRegion<2>& region);

/** Given a list of pixels, form them into an image, dilate the image,
  * and get the list of pixels in the dilated image. */
std::vector<itk::Index<2> > DilatePixelList(const std::vector<itk::Index<2> >& pixelList,
                                            const itk::ImageRegion<2>& region, const unsigned int radius);

/** Get the region of an image where a patch with radius 'patchRadius' centered on
  * every pixel will be entirely inside 'wholeRegion'. */
itk::ImageRegion<2> GetInternalRegion(const itk::ImageRegion<2>& wholeRegion, const unsigned int patchRadius);

/** Get all patches in 'imageRegion' with radius 'patchRadius' that contains a specified 'pixel'. */
std::vector<itk::ImageRegion<2> > GetAllPatchesContainingPixel(const itk::Index<2>& pixel,
                                                               const unsigned int patchRadius,
                                                               const itk::ImageRegion<2>& imageRegion);

/** Get all patches with radius 'patchRadius'. */
std::vector<itk::ImageRegion<2> > GetAllPatches(const itk::ImageRegion<2>& region, const unsigned int patchRadius);

/** Get the regions of patches surrounding every pixel in 'indices'. */
std::vector<itk::ImageRegion<2> > GetPatchesCenteredAtIndices(const std::vector<itk::Index<2> >& indices,
                                                              const unsigned int patchRadius);

/** Get the regions of patches surrounding every pixel in 'indices' if they are entirely inside 'imageRegion'. */
std::vector<itk::ImageRegion<2> > GetValidPatchesCenteredAtIndices(const std::vector<itk::Index<2> >& indices,
                                                                   const itk::ImageRegion<2>& imageRegion,
                                                                   const unsigned int patchRadius);

/** Find which location in 'pixels' is closest to 'queryPixel'. */
unsigned int ClosestIndexId(const std::vector<itk::Index<2> >& pixels, const itk::Index<2>& queryPixel);

/** Subtract 1 if necessary from each or either component to make both components even. */
itk::Size<2> MakeSizeEven(const itk::Size<2>& inputSize);

/** Find the distance between two indices. */
float IndexDistance(const itk::Index<2>& p0, const itk::Index<2>& p1);

/** Get all of the indices in a 'region'. */
std::vector<itk::Index<2> > GetIndicesInRegion(const itk::ImageRegion<2>& region);

/** Get a list of the indices in a 'region' downsampled in both dimensions by 'stride' . */
std::vector<itk::Index<2> > GetDownsampledIndicesInRegion(const itk::ImageRegion<2>& region,
                                                          const unsigned int stride);

/** Get the 4-neighbor indices around 'pixel' that are inside 'region'. */
std::vector<itk::Index<2> > Get4NeighborIndicesInsideRegion(const itk::Index<2>& pixel,
                                                            const itk::ImageRegion<2>& region);

/** Return the type of pixel that is in the file 'filename'. */
itk::ImageIOBase::IOComponentType GetPixelTypeFromFile(const std::string& filename);

/** Paraview requires 3D vectors to display glyphs, even if the vectors are really 2D.
    These functions appends a 0 to each vectors of a 2D vector image so that it can be
    easily visualized with Paraview. */
void Write2DVectorRegion(const FloatVector2ImageType* const image, const itk::ImageRegion<2>& region,
                         const std::string& filename);

/** Calls Write2DVectorRegion on a full image. */
void Write2DVectorImage(const FloatVector2ImageType* const image, const std::string& filename);

/** Determine if two pixels touch. */
bool IsNeighbor(const itk::Index<2>& index1, const itk::Index<2>& index2);

/** Dilate (expand) a region by a radius. */
itk::ImageRegion<2> DilateRegion(const itk::ImageRegion<2>& region, const unsigned int radius);

/** Erode (shrink) a region by a radius. */
itk::ImageRegion<2> ErodeRegion(const itk::ImageRegion<2>& region, const unsigned int radius);

/** Write an image where the pixels in 'regions' have been colored. */
void HighlightAndWriteRegions(const itk::Size<2>& imageSize, const std::vector<itk::ImageRegion<2> >& regions, const std::string& filename);

namespace detail
{
  /** Extract a channel of an image. The output image should be a scalar image,
  * but does not have to have the same pixel type as the input image. */
  template<typename TInputImage, typename TOutputImage>
  void ExtractChannel(const TInputImage* const image, const unsigned int channel,
                      TOutputImage* const output);
}

}// end namespace

#include "ITKHelpers.hpp"

#endif
