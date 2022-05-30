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
 * Scar Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkSurface.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkCleanPolyData.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>

// ITK
#include <itkPoint.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

// Qt
#include <QtDebug>
#include <QMessageBox>

// C++ Standard
#include <numeric>

// CemrgApp
#include "CemrgCommonUtils.h"
#include "CemrgScar3D.h"

CemrgScar3D::CemrgScar3D() {

    this->methodType = 2;
    this->minStep = -3, this->maxStep = 3;
    this->minScalar = 1E10, this->maxScalar = -1;
    this->voxelBasedProjection = false;
    this->debugging = false;
    this->scalars = vtkSmartPointer<vtkFloatArray>::New();
}

mitk::Surface::Pointer CemrgScar3D::Scar3D(std::string directory, mitk::Image::Pointer lgeImage, std::string segname) {

    //Convert to itk image
    itkImageType::Pointer scarImage;
    mitk::CastToItkImage(lgeImage, scarImage);
    itkImageType::Pointer visitedImage = itkImageType::New();
    ItkDeepCopy(scarImage, visitedImage);

    //Read in the mesh
    std::string path = directory + "/" + segname;
    mitk::Surface::Pointer surface = CemrgCommonUtils::LoadVTKMesh(path);
    vtkSmartPointer<vtkPo