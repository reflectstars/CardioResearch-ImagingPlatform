/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/*=========================================================================
CEMRG CMD APP TEMPLATE
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkSurface.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkCommandLineParser.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkImplicitBoolean.h>
#include <vtkImplicitVolume.h>
#include <vtkImplicitDataSet.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkDataSetSurfaceFilter.h>

// ITK
#include <itkPoint.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>

// Qt
#include <QtDebug>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>

// C++ Standard
#include <algorithm>
#include <string>
#include <numeric>

// CemrgApp
#include <CemrgScar3D.h>
#include <CemrgCommandLine.h>

typedef itk::Image<uint8_t, 3> ImageType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NNInterpolatorType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdType;
typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StrElType;
typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StrElType> ImFilterType;
typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> ImMultiplyType;
typedef itk::ImageRegionIterator<ImageType> IteratorType;

// ResampleImageFilterType::Pointer & imresize(ImageType::Pointer imToResize, ImageType::Pointer imTarget);
ThresholdType::Pointer thresholdImage(ImageType::Pointer input, uint8_t thresholdVal, QString debugPrefix, bool debugging);
ImFilterType::Pointer imdilate(ImageType::Pointer input, uint8_t radius, QString debugPrefix, bool debugging);
ImMultiplyType::Pointer immultiply(ImageType::Pointer input1, ImageType::Pointer input2, QString debugPrefix, bool debugging);
void relabel(ImageType::Pointer & imageToRelabel, ImageType::Pointer imageFragment, uint8_t fromLabel, uint8_t toLabel);

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("EASI processing");
    parser.setTitle("Labelled image re-labelling Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Relabel segmentation (RV/LV) for Fibres.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::InputFile,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
        "input", "i", mitkCommandLineParser::InputFile,
        "Input image segmentation", "Full path of segmentation file.",
        us::Any(), false);
    parser.addArgument(
        "bloodpool", "b", mitkCommandLineParser::InputFile,
        "Input image bloodpool seg", "Full path of bloodpool file.",
        us::Any(), false);
    parser.addArgument( // optional
        "bp-rv", "rv", mitkCommandLineParser::String,
        "Bloodpool segmentation label (RV)", "RV segmentation label (default=30)");
    parser.addArgument( // optional
        "bp-lv", "lv", mitkCommandLineParser::String,
        "Bloodpool segmentation label (LV)", "LV segmentation label (default=10)");
    parser.addArgument( // optional
        "swap-labels", "sl", mitkCommandLineParser::String,
        "Swap levels code", "Code to swap labels (syntax: 'N,M' , default:1,5)");
    parser.addArgument(// optional
        "output", "o", mitkCommandLineParser::String,
        "Output file name", "Where to save the output (default: output.nii).");
    parser.addArgument( // optional
        "debug", "d", mitkCommandLineParser::Bool,
        "Debug Outputs", "Save intermediate steps of processing.");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["input"].Empty() || parsedArgs["bloodpool"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);
    auto bloodFilename = us::any_cast<std::string>(parsedArgs["bloodpool"]);
    // auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    // Default values for optional arguments
    auto verbose = false;
    auto debug = false;
    std::string bp_rv = "30";
    std::string bp_lv = "10";
    std::string swaplabels = "1,5";
    std::string outFilename = "output.nii";

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }

    if (parsedArgs.end() != parsedArgs.find("bp-rv")) {
        bp_rv = us::any_cast<std::string>(parsedArgs["bp-rv"]);
    }

    if (parsedArgs.end() != parsedArgs.find("bp-lv")) {
        bp_lv = us::any_cast<std::string>(parsedArgs["bp-lv"]);
    }

    if (parsedArgs.end() != parsedArgs.find("swap-labels")) {
        swaplabels = us::any_cast<std::string>(parsedArgs["swap-labels"]);
    }

    if (parsedArgs.end() != parsedArgs.find("output")) {
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    }

    if (parsedArgs.end() != parsedArgs.find("debug")) {
        debug = us::any_cast<bool>(parsedArgs["debug"]);
    }


    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inname = QString::fromStdString(inFilename);
        QString bloodpname = QString::fromStdString(bloodFilename);
        QString outname = QString::fromStdStr