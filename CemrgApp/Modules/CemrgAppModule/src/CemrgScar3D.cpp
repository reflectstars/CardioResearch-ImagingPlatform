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
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

    //Calculate normals
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    vtkSmartPointer<vtkPolyData> tempPD = vtkSmartPointer<vtkPolyData>::New();
    tempPD->DeepCopy(pd);
    normals->ComputeCellNormalsOn();
    normals->SetInputData(tempPD);
    normals->SplittingOff();
    normals->Update();
    pd = normals->GetOutput();

    //Declarations
    vtkIdType numCellPoints;
    vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkFloatArray> cellNormals = vtkFloatArray::SafeDownCast(pd->GetCellData()->GetNormals());
    std::vector<double> allScalarsInShell;
    vtkSmartPointer<vtkFloatArray> scalarsOnlyStDev = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> scalarsOnlyMultiplier = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> scalarsOnlyIntensity = vtkSmartPointer<vtkFloatArray>::New();
    itkImageType::IndexType pixelXYZ;
    itkImageType::PointType pointXYZ;

    double maxSdev = -1e9;
    double maxSratio = -1e9;
    double mean = 0, var = 1;

    for (int i = 0; i < pd->GetNumberOfCells(); i++) {
        double pN[3];
        cellNormals->GetTuple(i, pN);
        double cX = 0, cY = 0, cZ = 0, numPoints = 0;
        pd->GetCellPoints(i, cellPoints);
        numCellPoints = cellPoints->GetNumberOfIds();

        for (vtkIdType neighborPoint = 0; neighborPoint < numCellPoints; ++neighborPoint) {

            //Get the neighbor point ID
            vtkIdType neighborPointID = cellPoints->GetId(neighborPoint);

            //Get the neighbor point position
            double cP[3];
            pd->GetPoint(neighborPointID, cP);

            // ITK method
            pointXYZ[0] = cP[0];
            pointXYZ[1] = cP[1];
            pointXYZ[2] = cP[2];
            scarImage->TransformPhysicalPointToIndex(pointXYZ, pixelXYZ);
            cP[0] = pixelXYZ[0];
            cP[1] = pixelXYZ[1];
            cP[2] = pixelXYZ[2];

            cX += cP[0];
            cY += cP[1];
            cZ += cP[2];

            numPoints++;
        }//_innerLoop

        cX /= numPoints;
        cY /= numPoints;
        cZ /= numPoints;

        // ITK method
        pointXYZ[0] = pN[0];
        pointXYZ[1] = pN[1];
        pointXYZ[2] = pN[2];
        scarImage->TransformPhysicalPointToIndex(pointXYZ, pixelXYZ);
        pN[0] = pixelXYZ[0];
        pN[1] = pixelXYZ[1];
        pN[2] = pixelXYZ[2];

        double scalar = GetIntensityAlongNormal(scarImage, visitedImage, pN[0], pN[1], pN[2], cX, cY, cZ);

        if (scalar > maxScalar) maxScalar = scalar;
        if (scalar < minScalar) minScalar = scalar;

        double sdev = (scalar - mean) / sqrt(var);
        // double sratio = mean ? scalar / mean : std::numeric_limits<decltype(sratio)>::max(); // mean=0 always false
        double sratio = std::numeric_limits<decltype(sratio)>::max(); // mean=0 always false

        if (maxSdev < sdev) maxSdev = sdev;
        if (maxSratio < sratio) maxSratio