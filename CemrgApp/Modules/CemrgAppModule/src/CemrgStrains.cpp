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
 * Strain Calculations Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkProperties.h>

// Qt
#include <QMessageBox>
#include <QDebug>

// Vtk
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkLineSource.h>
#include <vtkPlaneSource.h>
#include <vtkProbeFilter.h>
#include <vtkRegularPolygonSource.h>

// C++ Standard
#include <numeric>

// CemrgApp
#include "CemrgCommonUtils.h"
#include "CemrgStrains.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief TESTS remove later
 */
CemrgStrains::CemrgStrains() {
}

CemrgStrains::CemrgStrains(QString dir, int refMeshNo) {

    this->projectDirectory = dir;
    this->refAhaArea.assign(16, 0);
    this->refSurface = ReadVTKMesh(refMeshNo);
    this->refCellLabels.assign(refSurface->GetVtkPolyData()->GetNumberOfCells(), 0);
    this->refPointLabels.assign(refSurface->GetVtkPolyData()->GetNumberOfPoints(), 0.0);
    this->flatSurfScalars = vtkSmartPointer<vtkFloatArray>::New();
}

CemrgStrains::~CemrgStrains() {

    this->refArea.clear();
    this->refAhaArea.clear();
    this->refCellLabels.clear();
    this->refPointLabels.clear();
}

double CemrgStrains::CalculateGlobalSqzPlot(int meshNo) {

    //We want to load the mesh and then calculate the area
    mitk::Surface::Pointer refSurf = ReadVTKMesh(0);
    vtkSmartPointer<vtkPolyData> refPD = refSurf->GetVtkPolyData();
    mitk::Surface::Pointer surf = ReadVTKMesh(meshNo);
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();

    //Calculate squeeze
    double sqzValues = 0.0;
    for (vtkIdType cellID = 0; cellID < pd->GetNumberOfCells(); cellID++) {

        double refArea = GetCellArea(refPD, cellID);
        double area = GetCellArea(pd, cellID);
        double sqze = (area - refArea) / refArea;
        double wsqz = area * sqze;
        sqzValues += wsqz;

    }//_for

    //Average over entire mesh
    double avgSqzValues = sqzValues / pd->GetNumberOfCells();

    return avgSqzValues;
}

std::vector<double> CemrgStrains::CalculateSqzPlot(int meshNo) {

    if (refCellLabels.empty())
        return std::vector<double>(0);

    //We want to load the mesh and then calculate the area
    mitk::Surface::Pointer surf = ReadVTKMesh(meshNo);
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
    vtkSmartPointer<vtkFloatArray> sqzValues = vtkSmartPointer<vtkFloatArray>::New();

    //Calculate squeeze
    int index = 0;
    std::vector<double> squeeze(16, 0);
    for (vtkIdType cellID = 0; cellID < pd->GetNumberOfCells(); cellID++) {

        //Ignore non AHA segments
        if (refCellLabels[cellID] == 0) {
            sqzValues->InsertNextTuple1(0);
            continue;
        }//_if

        double area = GetCellArea(pd, cellID);
        double sqze = (area - refArea.at(index)) / refArea.at(index);
        double wsqz = area * sqze;
        squeeze.at(refCellLabels[cellID] - 1) += wsqz;

        //Global maps
        flatSurfScalars->InsertTuple1(index, wsqz);
        sqzValues->InsertNextTuple1(wsqz);

        index++;
    }//_for

    //Store squeeze values
    sqzValues->SetName("Squeez");
    pd->GetCellData()->AddArray(sqzValues);
    surf->SetVtkPolyData(pd);
    // QString path = projectDirectory + "/sqz-" + QString::number(meshNo) + ".vtk";
    // mitk::IOUtil::Save(surf, path.toStdString());

    //Average over AHA segments
    for (int i = 0; i < 16; i++)
        squeeze.at(i) /= refAhaArea.at(i);

    return squeeze;
}

std::vector<double> CemrgStrains::CalculateStrainsPlot(int meshNo, mitk::DataNode::Pointer lmNode, int flag) {

    if (refCellLabels.empty())
        return std::vector<double>(0);

    //We want to load the mesh and then calculate the strain
    std::vector<mitk::Point3D> lm = ConvertMPS(lmNode);
    mitk::Surface::Pointer surf = ReadVTKMesh(meshNo); ZeroVTKMesh(lm.at(0), surf);

    mitk::Point3D RIV2, centre;
    // Only do this for the manually marked landmark points (ap_3mv_2rv.mps)
    if (lm.size() == 6) {
        RIV2 = lm.at(5);
        centre = Circlefit3d(ZeroPoint(lm.at(0), lm.at(1)), ZeroPoint(lm.at(0), lm.at(2)), ZeroPoint(lm.at(0), lm.at(3)));
    } else {
        RIV2 = lm.at(3);
        centre = ZeroPoint(lm.at(0), lm.at(1));
    }

    mitk::Matrix<double, 3, 3> rotationMat = CalcRotationMatrix(centre, ZeroPoint(lm.at(0), RIV2)); RotateVTKMesh(rotationMat, surf);
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();

    //Radial, Circumferential, and Longitudinal strains for each AHA segment
    int index = 0;
    std::vector<double> strainRCL(16, 0);

    for (vtkIdType cellID = 0; cellID < pd->GetNumberOfCells(); cellID++) {

        //Ignore non AHA segments
        if (refCellLabels[cellID] == 0)
            continue;

        //Three nodes of the triangle
        vtkSmartPointer<vtkCell> cell = pd->GetCell(cellID);
        vtkSmartPointer<vtkTriangle> triangle = dynamic_cast<vtkTriangle*>(cell.GetPointer());
        double pt1[3], pt2[3], pt3[3];
        triangle->GetPoints()->GetPoint(0, pt1);
        triangle->GetPoints()->GetPoint(1, pt2);
        triangle->GetPoints()->GetPoint(2, pt3);

        //Coordinate system: vectors of the triangle
        mitk::Point3D vc1, vc2, vc3;
        vc1.SetElement(0, pt2[0] - pt1[0]);
        vc1.SetElement(1, pt2[1] - pt1[1]);
        vc1.SetElement(2, pt2[2] - pt1[2]);
        vc2.SetElement(0, pt3[0] - pt1[0]);
        vc2.SetElement(1, pt3[1] - pt1[1]);
        vc2.SetElement(2, pt3[2] - pt1[2]);
        vc3.SetElement(0, Cross(vc1, vc2).GetElement(0) / Norm(Cross(vc1, vc2)));
        vc3.SetElement(1, Cross(vc1, vc2).GetElement(1) / Norm(Cross(vc1, vc2)));
        vc3.SetElement(2, Cross(vc1, vc2).GetElement(2) / Norm(Cross(vc1, vc2)));

        //Assemble K matrix
        mitk::Matrix<double, 3, 3> K;
        K[0][0] = vc1.GetElement(0);
        K[0][1] = vc2.GetElement(0);
        K[0][2] = vc3.GetElement(0);
        K[1][0] = vc1.GetElement(1);
        K[1][1] = vc2.GetElement(1);
        K[1][2] = vc3.GetElement(1);
        K[2][0] = vc1.GetElement(2);
        K[2][1] = vc2.GetElement(2);
        K[2][2] = vc3.GetElement(2);

        //Calculate deformation gradient
        mitk::Matrix<double, 3, 3> F;
        F = K * refJ.at(index).GetInverse();

        //Calculate Strain Tensors
        mitk::Matrix<double, 3, 3> ET;
        mitk::Matrix<double, 3, 3> EYE; EYE.SetIdentity();

        //Green-Lagrange or Engineering
        if (flag > 2)
            ET = 0.5 * (F.GetTranspose() * F.GetVnlMatrix() - EYE.GetVnlMatrix());
        else
            ET = 0.5 * (F.GetVnlMatrix() + F.GetTranspose()) - EYE.GetVnlMatrix();

        //Rotate Green-Lagrange strain
        mitk::Matrix<double, 3, 3> E;
        E = refQ.at(index) * ET * refQ.at(index).GetTranspose();

        //Store symmetric Green-Lagrange strain in Voigt notation
        mitk::Matrix<double, 1, 3> EV;
        EV[0][0] = E(0, 0);
        EV[0][1] = E(1, 1);
        EV[0][2] = E(2, 2);

        //Prepare plot values
        strainRCL.at(refCellLabels[cellID] - 1) += EV[0][(flag > 2) ? flag - 2 : flag];
        flatSurfScalars->InsertTuple1(index, EV[0][(flag > 2) ? flag - 2 : flag]);

        index++;
    }//_for

    for (int i = 0; i < 16; i++)
        strainRCL.at(i) /= std::count(refCellLabels.begin(), refCellLabels.end(), i + 1);

    return strainRCL;
}

double CemrgStrains::CalculateSDI(std::vector<std::vector<double>> valueVectors, int cycleLengths, int noFrames) {

    if (valueVectors.size() == 0)
        return 0.0;

    std::vector<double> T2Ps(16, 0.0);
    std::vector<double> values(noFrames, 0.0);

    //Find time to peak values
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < noFrames; j++)
            values[j] = valueVectors[j][i];
        int index = std::distance(values.begin(), std::min_element(values.begin(), values.end()));
        double time2peak = index * (cycleLengths / noFrames);
        //T2Ps as a percentage of the cardiac cycle
        double t2p = (time2peak * 100) / cycleLengths;
        T2Ps.push_back(t2p);
    }

    //SD
    double sumDeviation = 0.0;
    double sum = std::accumulate(T2Ps.begin(), T2Ps.end(), 0.0);
    double mean = sum / T2Ps.size();
    for (unsigned int i = 0; i < T2Ps.size(); i++)
        sumDeviation += (T2Ps[i] - mean) * (T2Ps[i] - mean);
    return std::sqrt(sumDeviation / T2Ps.size());
}

std::vector<mitk::Surface::Pointer> CemrgStrains::ReferenceGuideLines(mitk::DataNode::Pointer lmNode) {

    //Prepare landmarks
    std::vector<mitk::Point3D> LandMarks = ConvertMPS(lmNode);
    mitk::Point3D CNTR, APEX, RIV1, RIV2, MIV1, MIV2, MIV3;

    if (LandMarks.size() >= 6) {

        APEX = LandMarks.at(0);
        MIV1 = LandMarks.at(1);
        MIV2 = LandMarks.at(2);
        MIV3 = LandMarks.at(3);
        RIV1 = LandMarks.at(4);
        RIV2 = LandMarks.at(5);
        //Calcaulte a circle through the mitral valve points
        CNTR = Circlefit3d(MIV1, MIV2, MIV3);

    } else if (LandMarks.size() == 4) {

        APEX = LandMarks.at(0);
        CNTR = LandMarks.at(1);
        RIV1 = LandMarks.at(2);
        RIV2 = LandMarks.at(3);
    }

    //Calculate points
    mitk::Point3D A1, A2, AB;
    mitk::Point3D PNT1, PNT2;
    A1 = RIV1 - APEX;
    A2 = RIV2 - APEX;
    AB = CNTR - APEX;
    PNT1.SetElement(0, APEX.GetElement(0) + Dot(A1, AB) / Dot(AB, AB) * AB.GetElement(0));
    PNT1.SetElement(1, APEX.GetElement(1) + Dot(A1, AB) / Dot(AB, AB) * AB.GetElement(1));
    PNT1.SetElement(2, APEX.GetElement(2) + Dot(A1, AB) / Dot(AB, AB) * AB.GetElement(2));
    PNT2.SetElement(0, APEX.GetElement(0) + Dot(A2, AB) / Dot(AB, AB) * AB.GetElement(0));
    PNT2.SetElement(1, APEX.GetElement(1) + Dot(A2, AB) / Dot(AB, AB) * AB.GetElement(1));
    PNT2.SetElement(2, APEX.GetElement(2) + Dot(A2, AB) / Dot(AB, AB) * AB.GetElement(2));

    //Draw guidelines
    mitk::Surface::Pointer line1 = mitk::Surface::New();
    vtkSmartPointer<vtkLineSource> lineSource1 = vtkSmartPointer<vtkLineSource>::New();
    lineSource1->SetPoint1(APEX.GetElement(0), APEX.GetElement(1), APEX.GetElement(2));
    lineSource1->SetPoint2(CNTR.GetElement(0), CNTR.GetElement(1), CNTR.GetElement(2));
    lineSource1->Update();
    line1->SetVtkPolyData(lineSource1->GetOutput());

    mitk::Surface::Pointer line2 = mitk::Surface::New();
    vtkSmartPointer<vtkLineSource> lineSource2 = vtkSmartPointer<vtkLineSource>::New();
    lineSource2->SetPoint1(RIV1.GetElement(0), RIV1.GetEleme