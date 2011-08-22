/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEntropyImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009-04-06 00:19:17 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEntropyImageFilter_txx
#define __itkEntropyImageFilter_txx

#include "itkEntropyImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
EntropyImageFilter<TInputImage, TOutputImage>
::EntropyImageFilter()
{
  m_Radius.Fill(1);
  m_NumberOfBins = 128;
  m_MarginalScale = 100;
  m_Mask = NULL;
}

template <class TInputImage, class TOutputImage>
void 
EntropyImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
  
  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_Radius );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    
    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

template< class TInputImage, class TOutputImage>
double
EntropyImageFilter< TInputImage, TOutputImage>
::CalculateEntropy(typename ListSampleType::Pointer listSample)
{
  // Calculate entropy from listSample  
  typename GeneratorType::Pointer generator = GeneratorType::New();
  typename GeneratorType::HistogramType::SizeType histSize;
  histSize.Fill(m_NumberOfBins); // number of bins
  generator->SetNumberOfBins(histSize);
  generator->SetAutoMinMax(true);
  generator->SetMarginalScale(m_MarginalScale);
  generator->SetListSample(listSample);
  generator->Update();

  typedef typename GeneratorType::HistogramType HistogramType;
  typename HistogramType::ConstPointer histogram = generator->GetOutput();
  float sum = histogram->GetTotalFrequency();
  float probability = 0;
  float entropy = 0;
  
  // Iterate through histogram and calculate entropy.
  typename HistogramType::ConstIterator itr = histogram->Begin();
  typename HistogramType::ConstIterator end = histogram->End();
  
  // Calculate entropy from histogram
  while(itr != end)
    {
    probability = itr.GetFrequency() / sum;
    if(probability > 0.99 / sum)
      {
      entropy += probability * log(probability) / log(2.0);
      }
    ++itr;
    }
   
  entropy = std::fabs(entropy);
  return entropy;
}

template< class TInputImage, class TOutputImage>
void
EntropyImageFilter< TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
  ConstantBoundaryCondition<InputImageType> nbc;
  nbc.SetConstant(0);
  
  ZeroFluxNeumannBoundaryCondition<InputImageType> zbc;

  ConstNeighborhoodIterator<InputImageType> bit;
  ImageRegionIterator<OutputImageType> it;
  
  // Allocate output
  typename OutputImageType::Pointer output = this->GetOutput();
  typename  InputImageType::ConstPointer input  = this->GetInput();
  
  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> bC;
  faceList = bC(input, outputRegionForThread, m_Radius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  
  double entropy;
  unsigned int num;
  
  // Process each of the boundary faces.  These are N-d regions which border
  // the edge of the buffer.
  for (fit=faceList.begin(); fit != faceList.end(); ++fit)
    { 
    if(m_Mask.IsNotNull())
      {
      bit = ConstNeighborhoodIterator<InputImageType>(m_Radius,
                                                      m_Mask, *fit);
      unsigned int neighborhoodSize = bit.Size();
      bit.OverrideBoundaryCondition(&nbc);
      bit.GoToBegin();                
      it = ImageRegionIterator<OutputImageType>(output, *fit);
      it.GotToBegin();
      while(!bit.IsAtEnd())
        {
        entropy = NumericTraits<double>::Zero;
        typename ListSampleType::Pointer listSample = ListSampleType::New();
        num = 0;
      
        if(bit.GetPixel(0) != 0)  // don't bother if neighborhood is at non-zero mask pixel
          {
          for(unsigned int i = 0; i < neighborhoodSize; ++i) // get number of non-zero voxels
            {
            if(bit.GetPixel(i) != 0)
              ++num;
            }
          listSample->SetMeasurementVectorSize(num); // set list sample
          for (unsigned int i = 0; i < neighborhoodSize; ++i)
            {
            if(bit.GetPixel(i) != 0) // only calculate entropy for non-zero mask pixels
              {
              listSample->PushBack(static_cast<InputPixelType>( input->GetPixel(bit.GetIndex(i)) ));
              }
            }
          
          // calculate the Entropy value
          entropy = CalculateEntropy(listSample);
          it.Set( static_cast<OutputPixelType>(entropy) );
          
          } // end if
        ++bit;
        ++it;
        progress.CompletedPixel();
        } // end while 
      } // end if mask
    else
      {
      bit = ConstNeighborhoodIterator<InputImageType>(m_Radius,
                                                    input, *fit);
      unsigned int neighborhoodSize = bit.Size();
      bit.OverrideBoundaryCondition(&zbc);
      bit.GoToBegin();                
      it = ImageRegionIterator<OutputImageType>(output, *fit);
      it.GotToBegin();
      while(!bit.IsAtEnd())
        {
        entropy = NumericTraits<double>::Zero;
        typename ListSampleType::Pointer listSample = ListSampleType::New();
        listSample->SetMeasurementVectorSize(neighborhoodSize);
        for (unsigned int i = 0; i < neighborhoodSize; ++i)
          {
          listSample->PushBack(static_cast<InputPixelType>( bit.GetPixel(i) ));
          }
       
        // calculate the Entropy value
        entropy = CalculateEntropy(listSample);
        it.Set( static_cast<OutputPixelType>(entropy) );
      
        ++bit;
        ++it;
        progress.CompletedPixel();
        }
      }
    }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput>
void
EntropyImageFilter<TInputImage, TOutput>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius: " << m_Radius << std::endl;
  os << indent << "Number of Bins: " << m_NumberOfBins << std::endl;
  os << indent << "Marginal Scale: " << m_MarginalScale << std::endl;

}

} // end namespace itk

#endif
