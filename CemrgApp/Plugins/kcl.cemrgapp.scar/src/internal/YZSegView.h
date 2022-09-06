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
 *
 * Scar Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * YZ
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef YZSegView_h
#define YZSegView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <mitkSurface.h>

// ITK
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImage.h>

#include <CemrgScar3D.h>
#include "ui_YZSegViewControls.h"

//Define the Image type and ITK iterator types
typedef itk::Image<short, 3> ImageType;
typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
typedef itk::ImageRegionIterator< ImageType> IteratorType;

/**
  \brief YZSegView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class YZSegView: public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;

    double FWMH_getPeak_SORT(ImageType::Pointer x);
    double GetStats(ImageType::Pointer im, ImageType::Pointer mask, int a);

    void thresholdImage(ImageType::Pointer img, ImageType::Pointer mask, double threshold);
    void insertion_sort(int array[], int l);


protected s