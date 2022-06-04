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
    lineSource2->SetPoint1(RIV1.GetElement(0), RIV1.GetElement(1), RIV1.GetElement(2));
    lineSource2->SetPoint2(PNT1.GetElement(0), PNT1.GetElement(1), PNT1.GetElement(2));
    lineSource2->Update();
    line2->SetVtkPolyData(lineSource2->GetOutput());

    mitk::Surface::Pointer line3 = mitk::Surface::New();
    vtkSmartPointer<vtkLineSource> lineSource3 = vtkSmartPointer<vtkLineSource>::New();
    lineSource3->SetPoint1(RIV2.GetElement(0), RIV2.GetElement(1), RIV2.GetElement(2));
    lineSource3->SetPoint2(PNT2.GetElement(0), PNT2.GetElement(1), PNT2.GetElement(2));
    lineSource3->Update();
    line3->SetVtkPolyData(lineSource3->GetOutput());

    /* // Commenting out the circle guideline as we dont have 3 MV points with Siemen Landmarks
    mitk::Surface::Pointer circle = mitk::Surface::New();
    vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
    plane->SetOrigin(CNTR.GetElement(0), CNTR.GetElement(1), CNTR.GetElement(2));
    plane->SetPoint1(MIV1.GetElement(0), MIV1.GetElement(1), MIV1.GetElement(2));
    plane->SetPoint2(MIV2.GetElement(0), MIV2.GetElement(1), MIV2.GetElement(2));
    plane->Update();
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
    polygonSource->GeneratePolygonOff();
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(CNTR.GetElement(0), CNTR.GetElement(1), CNTR.GetElement(2));
    polygonSource->SetNormal(plane->GetNormal());
    polygonSource->Update();
    circle->SetVtkPolyData(polygonSource->GetOutput()); //*/

    //Prepare guidelines vector
    std::vector<mitk::Surface::Pointer> guidelines;
    guidelines.push_back(line1);
    guidelines.push_back(line2);
    guidelines.push_back(line3);
    //guidelines.push_back(circle);
    return guidelines;
}

mitk::Surface::Pointer CemrgStrains::ReferenceAHA(mitk::DataNode::Pointer lmNode, int segRatios[], bool pacingSite) {

    //Read the mesh data
    vtkSmartPointer<vtkPolyData> pd = refSurface->GetVtkPolyData();

    //Prepare landmarks
    std::vector<mitk::Point3D> LandMarks = ConvertMPS(lmNode);
    // Siemens 4 LM = [apex, basecenter, RV1, RV2];
    // Siemens 7 LM = [apex, baseMV1, baseMV2, baseMV3, RV1, RV2, apex
    // (just to make it 7 and different from manual LM)]
    // Manual 6 LM = [apex, MV1, MV2, MV3, RV1, RV2]
    if (LandMarks.size() != 4 && LandMarks.size() != 6 && LandMarks.size() != 7)
        return refSurface;

    mitk::Point3D centre, RIV1, RIV2, APEX, MIV1, MIV2, MIV3;
    if (LandMarks.size() >= 6) {

        APEX = LandMarks.at(0);
        MIV1 = LandMarks.at(1);
        MIV2 = LandMarks.at(2);
        MIV3 = LandMarks.at(3);
        RIV1 = LandMarks.at(4);
        RIV2 = LandMarks.at(5);
        // Calcaulte a circle through the mitral valve points
        centre = Circlefit3d(ZeroPoint(APEX, MIV1), ZeroPoint(APEX, MIV2), ZeroPoint(APEX, MIV3));

    } else if (LandMarks.size() == 4) {

        // Zero all points relative to apex
        APEX = LandMarks.at(0);
        RIV1 = LandMarks.at(2);
        RIV2 = LandMarks.at(3);
        centre = ZeroPoint(APEX, LandMarks.at(1));
    }

    //Zero all points relative to apex
    RIV1 = ZeroPoint(APEX, RIV1);
    RIV2 = ZeroPoint(APEX, RIV2);
    ZeroVTKMesh(APEX, refSurface);
    APEX = ZeroPoint(APEX, APEX);

    //Calculate a circle through the mitral valve points
    mitk::Point3D RCTR;

    //Define Rotation matrix
    mitk::Matrix<double, 3, 3> rotationMat = CalcRotationMatrix(centre, RIV2);

    //Rotate points to new frame
    if (LandMarks.size() >= 6) {
        MIV1 = RotatePoint(rotationMat, MIV1);
        MIV2 = RotatePoint(rotationMat, MIV2);
        MIV3 = RotatePoint(rotationMat, MIV3);
    }
    RIV1 = RotatePoint(rotationMat, RIV1);
    RIV2 = RotatePoint(rotationMat, RIV2);
    RCTR = RotatePoint(rotationMat, centre);

    /**
      TEST
      **/
      //qDebug() << "RCTR IS " << RCTR.GetElement(0) << RCTR.GetElement(1) << RCTR.GetElement(2);

      //Rotate mesh to new frame
    RotateVTKMesh(rotationMat, refSurface);

    //Find the mesh Z range
    //double min = GetMinMax(pd,2).at(0);
    //double max = GetMinMax(pd,2).at(1);
    double min = APEX.GetElement(2);
    double max = RCTR.GetElement(2);
    double RangeZ = (max - min) * 1.0; //0.99;

    //Top, mid, and base segments heights
    double TOP = RangeZ * (segRatios[0] / 100.00 + segRatios[1] / 100.00 + segRatios[2] / 100.0) + min;
    double MID = RangeZ * (segRatios[1] / 100.00 + segRatios[2] / 100.0) + min;
    double BAS = RangeZ * (segRatios[2] / 100.0) + min;

    //Angle RV cusp 2
    double RVangle1 = atan2(RIV1.GetElement(1), RIV1.GetElement(0));
    double RVangle2 = atan2(RIV2.GetElement(1), RIV2.GetElement(0));

    //Assuming RV angle1 < RV angle 2
    double appendAngle;
    double sepA, freeA;

    if ((LandMarks.size() == 6) && (pacingSite == false)) {

        // only do this for manual segmentation, with 6 points
        sepA = (RVangle2 - RVangle1) / 2;
        freeA = (2 * M_PI - (RVangle2 - RVangle1)) / 4;
        appendAngle = -RVangle1;
        //appendAngle -= freeA - ( M_PI / 3 );

    } else {

        sepA = (2 * M_PI) / 6;
        freeA = (2 * M_PI) / 6;
        if (RVangle1 > 0) {
            appendAngle = -((M_PI - RVangle1) / 2 + RVangle1) + M_PI / 3;
        } else {
            appendAngle = -RVangle1 / 2 + M_PI / 3;
        }//_if
        if (pacingSite == true) {
            appendAngle += M_PI / 6;
        }

    }//_if
    qDebug() << "appendAngle " << appendAngle;

    //Point angles
    std::vector<double> pAngles;
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double pAngle;
        double* pt = pd->GetPoint(i);
        pAngle = atan2(pt[1], pt[0]) + appendAngle;
        pAngle = pAngle * (pAngle > 0 ? 1 : 0) + (2 * M_PI + pAngle) * (pAngle < 0 ? 1 : 0);
        pAngles.push_back(pAngle);
    }

    //Centre angles
    std::vector<double> cAngles;
    for (vtkIdType cellID = 0; cellID < pd->GetNumberOfCells(); cellID++) {
        double cAngle;
        mitk::Point3D ctrT = GetCellCenter(pd, cellID);
        cAngle = atan2(ctrT.GetElement(1), ctrT.GetElement(0)) + appendAngle;
        cAngle = cAngle * (cAngle > 0 ? 1 : 0) + (2 * M_PI + cAngle) * (cAngle < 0 ? 1 : 0);
        cAngles.push_back(cAngle);
    }

    //Base points
    std::vector<int> pBindex;
    std::vector<int> pMindex;
    std::vector<int> pAindex;
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double* pt = pd->GetPoint(i);
        pBindex.push_back((pt[2] >= MID ? 1 : 0) * (pt[2] <= TOP ? 1 : 0));
        pMindex.push_back((pt[2] >= BAS ? 1 : 0) * (pt[2] < MID ? 1 : 0));
        pAindex.push_back((pt[2] < BAS ? 1 : 0));
    }

    //Centre points
    std::vector<int> cBindex;
    std::vector<int> cMindex;
    std::vector<int> cAindex;
    for (vtkIdType cellID = 0; cellID < pd->GetNumberOfCells(); cellID++) {
        mitk::Point3D ctrT = GetCellCenter(pd, cellID);
        cBindex.push_back((ctrT.GetElement(2) >= MID ? 1 : 0) * (ctrT.GetElement(2) < TOP ? 1 : 0));
        cMindex.push_back((ctrT.GetElement(2) >= BAS ? 1 : 0) * (ctrT.GetElement(2) < MID ? 1 : 0));
        cAindex.push_back((ctrT.GetElement(2) < BAS ? 1 : 0));
    }

    //Assign points labels for 3 layers
    AssignpLabels(0, refPointLabels, pBindex, pAngles, sepA, freeA);
    AssignpLabels(1, refPointLabels, pMindex, pAngles, sepA, freeA);
    AssignpLabels(2, refPointLabels, pAindex, pAngles, sepA, freeA);
    AssigncLabels(0, refCellLabels, cBindex, cAngles, sepA, freeA);
    AssigncLabels(1, refCellLabels, cMindex, cAngles, sepA, freeA);
    AssigncLabels(2, refCellLabels, cAindex, cAngles, sepA, freeA);

    //Calculate reference mesh attributes
    for (vtkIdType cellID = 0; cellID < pd->GetNumberOfCells(); cellID++) {

        //Ignore non AHA segments
        if (refCellLabels[cellID] == 0)
            continue;

        //Area
        double area = GetCellArea(pd, cellID);
        refArea.push_back(area);
        refAhaArea.at(refCellLabels[cellID] - 1) += area;

        //Axis
        vtkSmartPointer<vtkCell> cell = pd->GetCell(cellID);
        mitk::Matrix<double, 3, 3> J;
        mitk::Matrix<double, 3, 3> Q = GetCellAxes(cell, RCTR, J);
        refJ.push_back(J);
        refQ.push_back(Q);
    }

    //Setup flattened AHA mesh
    flatSurface = refSurface->Clone();
    vtkSmartPointer<vtkPolyData> poly = flatSurface->GetVtkPolyData();
    for (int i = 0; i < poly->GetNumberOfPoints(); i++) {
        double* point = poly->GetPoint(i);
        double  radii = point[2];
        double  theta = pAngles.at(i);
        point[0] = radii * cos(theta);
        point[1] = radii * sin(theta);
        point[2] = 0;
        poly->GetPoints()->SetPoint(i, point);
    }//_for
    poly->BuildLinks();
    for (vtkIdType cellID = 0; cellID < poly->GetNumberOfCells(); cellID++)
        if (refCellLabels[cellID] == 0)
            poly->DeleteCell(cellID);
    poly->RemoveDeletedCells();

    //Setup colors of the AHA segmentations
    vtkSmartPointer<vtkUnsignedCharArray> segmentColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    segmentColors->SetNumberOfComponents(3);
    segmentColors->SetNumberOfTuples(pd->GetNumberOfPoints());
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        float rgbA[3];
        std::vector<float> rgbV = GetAHAColour((int)refPointLabels.at(i));
        std::copy(rgbV.begin(), rgbV.end(), rgbA);
        segmentColors->InsertTuple(i, rgbA);
    }
    refSurface->GetVtkPolyData()->GetPointData()->SetScalars(segmentColors);

    //Return to the original position and rotation
    APEX = ZeroPoint(LandMarks.at(0), APEX);
    RotateVTKMesh(rotationMat.GetInverse().as_matrix(), refSurface);
    ZeroVTKMesh(APEX, refSurface);
    return refSurface;
}

mitk::Surface::Pointer CemrgStrains::FlattenedAHA() {

    if (refCellLabels.empty())
        return mitk::Surface::New();
    return flatSurface;
}

vtkSmartPointer<vtkFloatArray> CemrgStrains::GetFlatSurfScalars() const {

    return flatSurfScalars;
}

std::vector<float> CemrgStrains::GetAHAColour(int label) {

    switch (label) {
    case 1:
        return std::vector<float>{230, 25, 75};
    case 2:
        return std::vector<float>{60, 180, 75};
    case 3:
        return std::vector<float>{255, 225, 25};
    case 4:
        return std::vector<float>{0, 130, 200};
    case 5:
        return std::vector<float>{245, 130, 48};
    case 6:
        return std::vector<float>{145, 30, 180};
    case 7:
        return std::vector<float>{70, 240, 240};
    case 8:
        return std::vector<float>{250, 190, 190};
    case 9:
        return std::vector<float>{0, 128, 128};
    case 10:
        return std::vector<float>{170, 110, 40};
    case 11:
        return std::vector<float>{128, 0, 0};
    case 12:
        return std::vector<float>{128, 128, 0};
    case 13:
        return std::vector<float>{0, 0, 128};
    case 14:
        return std::vector<float>{128, 128, 128};
    case 15:
        return std::vector<float>{255, 255, 255};
    case 16:
        return std::vector<float>{230, 190, 255};
    }
    return std::vector<float>{0, 0, 0};
}

/**************************************************************************************************
 *************** HELPER FUNCTIONS *****************************************************************
 **************************************************************************************************/

mitk::Point3D CemrgStrains::ZeroPoint(mitk::Point3D apex, mitk::Point3D point) {

    //Zero relative to the apex
    point.SetElement(0, point.GetElement(0) - apex.GetElement(0));
    point.SetElement(1, point.GetElement(1) - apex.GetElement(1));
    point.SetElement(2, point.GetElement(2) - apex.GetElement(2));
    return point;
}

void CemrgStrains::ZeroVTKMesh(mitk::Point3D apex, mitk::Surface::Pointer surface) {

    //Retrieve the data
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = point[0] - apex.GetElement(0);
        point[1] = point[1] - apex.GetElement(1);
        point[2] = point[2] - apex.GetElement(2);
        pd->GetPoints()->SetPoint(i, point);
    }
}

mitk::Point3D CemrgStrains::RotatePoint(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Point3D point) {

    mitk::Matrix<double, 1, 3> vec;
    mitk::Matrix<double, 3, 1> ans;

    vec[0][0] = point.GetElement(0);
    vec[0][1] = point.GetElement(1);
    vec[0][2] = point.GetElement(2);
    ans = rotationMatrix * vec.GetTranspose();
    point.SetElement(0, ans[0][0]);
    point.SetElement(1, ans[1][0]);
    point.SetElement(2, ans[2][0]);

    return point;
}

void CemrgStrains::RotateVTKMesh(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Surface::Pointer surface) {

    //Retrieve the data
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

    //Rotate mesh into new coordiantes
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        mitk::Point3D point;
        double* pt = pd->GetPoint(i);
        point.SetElement(0, pt[0]);
        point.SetElement(1, pt[1]);
        point.SetElement(2, pt[2]);
        point = RotatePoint(rotationMatrix, point);
        pt[0] = point.GetElement(0);
        pt[1] = point.GetElement(1);
        pt[2] = point.GetElement(2);
        pd->GetPoints()->SetPoint(i, pt);
    }
}

double CemrgStrains::GetCellArea(vtkSmartPointer<vtkPolyData> pd, vtkIdType cellID) {

    vtkSmartPointer<vtkCell> cell = pd->GetCell(cellID);
    vtkSmartPointer<vtkTriangle> triangle = dynamic_cast<vtkTriangle*>(cell.GetPointer());
    double pt1[3], pt2[3], pt3[3];
    triangle->GetPoints()->GetPoint(0, pt1);
    triangle->GetPoints()->GetPoint(1, pt2);
    triangle->GetPoints()->GetPoint(2, pt3);

    double area = vtkTriangle::TriangleArea(pt1, pt2, pt3);
    return area;
}

mitk::Point3D CemrgStrains::GetCellCenter(vtkSmartPointer<vtkPolyData> pd, vtkIdType cellID) {

    mitk::Point3D centre;
    vtkSmartPointer<vtkCell> cell = pd->GetCell(cellID);
    vtkSmartPointer<vtkTriangle> triangle = dynamic_cast<vtkTriangle*>(cell.GetPointer());
    double pt1[3], pt2[3], pt3[3];
    triangle->GetPoints()->GetPoint(0, pt1);
    triangle->GetPoints()->GetPoint(1, pt2);
    triangle->GetPoints()->GetPoint(2, pt3);

    double x = pt1[0] + pt2[0] + pt3[0];
    double y = pt1[1] + pt2[1] + pt3[1];
    double z = pt1[2] + pt2[2] + pt3[2];
    centre.SetElement(0, x / 3);
    centre.SetElement(1, y / 3);
    centre.SetElement(2, z / 3);
    return centre;
}

mitk::Matrix<double, 3, 3> CemrgStrains::GetCellAxes(vtkSmartPointer<vtkCell>& cell, mitk::Point3D& termPt, mitk::Matrix<double, 3, 3>& J) {

    //Three nodes of the triangle
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
    vc3