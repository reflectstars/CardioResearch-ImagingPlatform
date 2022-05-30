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
        if (maxSratio < sratio) maxSratio = sratio;

        /**
         * Placeholder variables for potential GUI options.
         * We removed them in favour of a simplified functionality
         *
        int _ONLY_POSITIVE_STDEVS = 1;
        int _SCAR_AS_STANDARD_DEVIATION = 1;
        int _SCAR_MIP = 1;
        *
        * Comments that refer to this are identified with the following:
        * [placeholder]
         *
         */

        // [placeholder] simplified functionality without placeholder variables
        if (sdev < 0) sdev = 0;

        // [placeholder] comlete functionality with placeholder variables
        // if (_ONLY_POSITIVE_STDEVS == 1 && sdev < 0) sdev = 0;

        scalarsOnlyStDev->InsertTuple1(i, sdev);
        scalarsOnlyIntensity->InsertTuple1(i, scalar);
        scalarsOnlyMultiplier->InsertTuple1(i, sratio);

        //For default scalar to plot
        double scalarToPlot = (scalar - mean) / sqrt(var);

        if (scalarToPlot <= 0) scalarToPlot = 0;

        // [placeholder] simplified functionality without placeholder variables
        scalars->InsertTuple1(i, scalarToPlot);
        allScalarsInShell.push_back(scalarToPlot);

        // [placeholder] comlete functionality with placeholder variables
        // if (_SCAR_MIP == 1 && _SCAR_AS_STANDARD_DEVIATION == 1) {
        //      scalars->InsertTuple1(i, scalarToPlot);
        //      allScalarsInShell.push_back(scalarToPlot);
        // } else {
        //      scalars->InsertTuple1(i, scalar);
        //      allScalarsInShell.push_back(scalar); }
    }//_for

    scarDebugLabel = visitedImage;
    pd->GetCellData()->SetScalars(scalars);
    surface->SetVtkPolyData(pd);
    return surface;
}

mitk::Surface::Pointer CemrgScar3D::ClipMesh3D(mitk::Surface::Pointer surface, mitk::PointSet::Pointer landmarks) {

    //Retrieve mean and distance of 3 points
    double x_c = 0;
    double y_c = 0;
    double z_c = 0;
    for (int i = 0; i < landmarks->GetSize(); i++) {
        x_c = x_c + landmarks->GetPoint(i).GetElement(0);
        y_c = y_c + landmarks->GetPoint(i).GetElement(1);
        z_c = z_c + landmarks->GetPoint(i).GetElement(2);
    }//_for
    x_c /= landmarks->GetSize();
    y_c /= landmarks->GetSize();
    z_c /= landmarks->GetSize();
    double * distance = new double[landmarks->GetSize()];
    for (int i = 0; i < landmarks->GetSize(); i++) {
        double x_d = landmarks->GetPoint(i).GetElement(0) - x_c;
        double y_d = landmarks->GetPoint(i).GetElement(1) - y_c;
        double z_d = landmarks->GetPoint(i).GetElement(2) - z_c;
        distance[i] = sqrt(pow(x_d, 2) + pow(y_d, 2) + pow(z_d, 2));
    }//_for
    double radius = *std::max_element(distance, distance + landmarks->GetSize());
    double centre[3] = {x_c, y_c, z_c};

    //Clipper
    vtkSmartPointer<vtkSphere> sphere = vtkSmartPointer<vtkSphere>::New();
    sphere->SetCenter(centre);
    sphere->SetRadius(radius);
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetClipFunction(sphere);
    clipper->SetInputData(surface->GetVtkPolyData());
    clipper->InsideOutOff();
    clipper->Update();

    //Extract and clean surface mesh
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfer->SetInputData(clipper->GetOutput());
    surfer->Update();
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(surfer->GetOutputPort());
    cleaner->Update();
    vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    lrgRegion->SetInputConnection(cleaner->GetOutputPort());
    lrgRegion->SetExtractionModeToLargestRegion();
    lrgRegion->Update();
    cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(lrgRegion->GetOutputPort());
    cleaner->Update();

    //Return the clipped mesh
    surface->SetVtkPolyData(cleaner->GetOutput());
    return surface;
}

bool CemrgScar3D::CalculateMeanStd(mitk::Image::Pointer lgeImage, mitk::Image::Pointer roiImage, double& mean, double& stdv) {

    //Access image volumes
    mitk::ImagePixelReadAccessor<float, 3> readAccess1(lgeImage);
    float* pvLGE = (float*)readAccess1.GetData();
    mitk::ImagePixelReadAccessor<float, 3> readAccess2(roiImage);
    float* pvROI = (float*)readAccess2.GetData();

    int dimsLGE = lgeImage->GetDimensions()[0] * lgeImage->GetDimensions()[1] * lgeImage->GetDimensions()[2];
    int dimsROI = roiImage->GetDimensions()[0] * roiImage->GetDimensions()[1] * roiImage->GetDimensions()[2];
    if (dimsLGE != dimsROI) {
        QMessageBox::critical(NULL, "Attention", "The mask and the image dimensions do not match!");
        return false;
    }//_wrong dimensions

    //Loop image voxels
    std::vector<float> voxelValues;
    for (int i = 0; i < dimsROI; i++) {
        if (*pvROI == 1)
            voxelValues.push_back(*pvLGE);
        pvLGE++;
        pvROI++;
    }//_for

    //Calculate mean and std
    double sumDeviation = 0.0;
    double sum = std::accumulate(voxelValues.begin(), voxelValues.end(), 0.0);
    mean = sum / voxelValues.size();
    for (unsigned int i = 0; i < voxelValues.size(); i++)
        sumDeviation += (voxelValues[i] - mean) * (voxelValues[i] - mean);
    stdv = std::sqrt(sumDeviation / voxelValues.size());
    return true;
}

double CemrgScar3D::Thresholding(double thresh) {

    int ctr1 = 0, ctr2 = 0;
    for (int i = 0; i < scalars->GetNumberOfTuples(); i++) {
        double value = scalars->GetValue(i);
        if (value == -1) {
            ctr1++;
            continue;
        }//_if
        if (value > thresh) ctr2++;
    }
    double percentage = (ctr2 * 100.0) / (scalars->GetNumberOfTuples() - ctr1);
    return percentage;
}

void CemrgScar3D::SaveNormalisedScalars(double divisor, mitk::Surface::Pointer surface, QString name) {

    MITK_INFO << "Dividing by the mean value of the bloodpool.";
    MITK_ERROR(divisor == 0) << "Can't divide by 0.";

    vtkSmartPointer<vtkCellDataToPointData> cell_to_point = vtkSmartPointer<vtkCellDataToPointData>::New();
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
    vtkSmartPointer<vtkFloatArray> normalisedScalars = vtkSmartPo