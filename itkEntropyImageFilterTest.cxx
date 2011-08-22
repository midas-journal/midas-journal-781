/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEntropyImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007-08-10 14:34:02 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <fstream>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPNGImageIOFactory.h"
#include "itkTextOutput.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNumericTraits.h"
#include "itkFilterWatcher.h"
#include "itkEntropyImageFilter.h"

int main(int ac, char* av[] )
{
  // Comment the following if you want to use the itk text output window
  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  if(ac < 6)
    {
    std::cerr << "Usage: " << av[0] << " InputImage InputImageExtension Radius NumberOfBins MarginalScale\n";
    return -1;
    }

  itk::Size<2> radius;
  typedef itk::Image<unsigned char, 2> myImageIn;
  typedef itk::Image<float, 2> myImageOut;
  typedef itk::Image<unsigned char, 2> myImageChar;
  itk::ImageFileReader<myImageIn>::Pointer input = itk::ImageFileReader<myImageIn>::New();
  input->SetFileName(std::string(av[1]) + std::string(av[2]));
  
  // Create a filter
  typedef itk::EntropyImageFilter<myImageIn, myImageOut> FilterType;
  FilterType::Pointer filter = FilterType::New();
  FilterWatcher filterWatch(filter);

  typedef itk::RescaleIntensityImageFilter<myImageOut, myImageChar> RescaleFilterType; 
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetOutputMinimum(0);
  rescale->SetOutputMaximum(255);
  rescale->SetInput(filter->GetOutput());

  try
    {
    radius.Fill(std::atoi(av[3]));
    filter->SetInput(input->GetOutput());
    filter->SetRadius(radius);
    filter->SetNumberOfBins(std::atoi(av[4]));
    filter->SetMarginalScale(std::atoi(av[5]));
    filter->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception detected: "  << e.GetDescription();
    return -1;
    }

  // Write output as float
  itk::ImageFileWriter<myImageOut>::Pointer writer;
  writer = itk::ImageFileWriter<myImageOut>::New();
  writer->SetInput(filter->GetOutput());
  writer->SetFileName(std::string(av[1]) + std::string("_out.nii"));
  writer->Update();  
    
  // Write output cast to an image of unsigned char
  itk::ImageFileWriter<myImageChar>::Pointer writer2;
  writer2 = itk::ImageFileWriter<myImageChar>::New();
  writer2->SetInput(rescale->GetOutput());
  writer2->SetFileName(std::string(av[1]) + std::string("_out.png"));
  writer2->Update();

  return EXIT_SUCCESS;
}
