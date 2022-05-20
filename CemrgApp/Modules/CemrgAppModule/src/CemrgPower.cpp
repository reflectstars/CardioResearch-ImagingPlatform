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
 * Power Transmitter Calculations Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * angela.lee@kcl.ac.uk
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
#include <QFileInfo>
#include <QCoreApplication>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkMatrix3x3.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkExtractSelection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkInformation.h>

// C++ Standard
#include <string.h>
#include <vector>


#include "CemrgPower.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

CemrgPower::CemrgPower() {
}

CemrgPower::CemrgPower(QString dir, int ribSpacing) {

    this->projectDirectory = dir;
    this->currentRibSpacing = ribSpacing;
}

mitk::Surface::Pointer CemrgPower::MapPowerTransmitterToLandmarks(mitk::DataNode::Pointer lmNode) {

    mitk::Surface::Pointer mesh;

    // Read EBR vtk mesh
    QString in_ebr_fname;
    in_ebr_fname = QCoreApplication::applicationDirPath() + "/EBR_data/ebr_initial.vtk";
    if (!QFileInfo::exists(in_ebr_fname)) {
        QMessageBox::warning(NULL, "Attention", "Power transmitter template does not exist: Please check that programPath/EBR_data/ebr_initial.vtk exists!");
        return mesh;
    }

    //QByteArray ba = in_ebr_fname.toLocal8Bit();
    //const char *c_str2 = ba.data();
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();

    reader->SetFileName(in_ebr_fname.toLocal8Bit().data());
    reader->Update();
    vtkPolyData* polydata = reader->GetOutput();

    // Read the location of the landmarks indicating the possible location of the power transmitter. In MITK plugin this will be a mps file, but in this separate executable can just have it as a points file
    //Prepare landmark
    std::vector<mitk::Point3D> mps = CemrgPower::ConvertMPS(lmNode);
    double vl_cas[3];
    double vl_ebr[3] = {0.985814, 0.107894, -0.128566};
    double vw_ebr[3] = {0.0193834, -0.277855, -0.960427};
    double x0[3] = {255.022, 60.4436, 230.751};

    /*
    // Hard coded indices on the EBR mesh that represent:
    polydata->GetPoint(113747,x0); // Mid medial edge of the transmitter (edge facing rib sternum)
    polydata->GetPoint(69982,x1); // Along the length
    polydata->GetPoint(112413,w1); // Along the width

    vl_ebr=norm(x1-x0);
    vw_ebr=norm(w1-x0);
    */

    // Check that there are more than 2 LM points selected
    if (mps.size() >= 2) {
        // For the MPS landmarks
        vl_cas[0] = (mps[1][0] - mps[0][0]);
        vl_cas[1] = (mps[1][1] - mps[0][1]);
        vl_cas[2] = (mps[1][2] - mps[0][2]);
        CemrgPower::normalise(vl_cas);
    } else {
        QMessageBox::warning(NULL, "Attention", "Select 2 or 3 landmark points");
    }
    // Find the rotation matrix for the vectors along length of transmitter (EBR template and MPS landmarks)
    //Matrix3f R=Matrix3f::SetIdentity();
    vtkSmartPointer<vtkMatrix3x3> R = vtkSmartPointer<vtkMatrix3x3>::New();
    // Get rotation matrix between template EBR transmitter and length vector from LMs
    CemrgPower::fcn_RotationFromTwoVectors(vl_ebr, vl_cas, R);

    vtkSmartPointer<vtkMatrix4x4> R1 = vtkSmartPointer<vtkMatrix4x4>::New();

    R1->SetElement(0, 0, R->GetElement(0, 0));
    R1->SetElement(0, 1, R->GetElement(0, 1));
    R1->SetElement(0, 2, R->GetElement(0, 2));
    R1->SetElement(1, 0, R->GetElement(1, 0));
    R1->SetElement(1, 1, R->GetElement(1, 1));
    R1->SetElement(1, 2, R->GetElement(1, 2));
    R1->SetElement(2, 0, R->GetElement(2, 0));
    R1->SetElement(2, 1, R->GetElement(2, 1));
    R1->SetElement(2, 2, R->GetElement(2, 2));

    double vw_ebr_rot[3];
    // Apply rotation matrix to vw_ebr vector
    R->MultiplyPoint(vw_ebr, vw_ebr_rot);

    vtkSmartPointer<vtkMatrix4x4> R2 = vtkSmartPointer<vtkMatrix4x4>::New();
    // Find the rotation matrix for the vectors along width of transmitter (EBR template and MPS landmarks)
    // This is an optional step if 3 landmark points are selected
    if (mps.size() >= 3) {
        // Re-initialise R;
        R = vtkSmartPointer<vtkMatrix3x3>::New();
        double vw_cas[3] = {mps[2][0] - mps[0][0], mps[2][1] - mps[0][1], mps[2][2] - mps[0][2]};
        CemrgPower::normalise(vw_cas);

        CemrgPower::fcn_RotationFromTwoVectors(vw_ebr_rot, vw_cas, R);

        R2->SetElement(0, 0, R->GetElement(0, 0));
        R2->SetElement(0, 1, R->GetElement(0, 1));
        R2->SetElement(0, 2, R->GetElement(0, 2));
        R2->SetElement(1, 0, R->GetElement(1, 0));
        R2->SetElement(1, 1, R->GetElement(1, 1));
        R2->SetElement(1, 2, R->GetElement(1, 2));
        R2->SetElement(2, 0, R->GetElement(2, 0));
        R2->SetElement(2, 1, R->GetElement(2, 1));
        R2->SetElement(2, 2, R->GetElement(2, 2));
    }
    // We want to apply a rotation around the location of pinned location in the rib space

    vtkSmartPointer<vtkMatrix4x4> vtk_T = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> T0 = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> T1 = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> T2 = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> T3 = vtkSmartPointer<vtkMatrix4x4>::New();

    T0->SetElement(0, 3, -x0[0]);
    T0->SetElement(1, 3, -x0[1]);
    T0->SetElement(2, 3, -x0[2]);

    vtkMatrix4x4::Multiply4x4(R2, T0, T1);
    vtkMatrix4x4::Multiply4x4(R1, T1, T2);

    T0->SetElement(0, 3, x0[0]);
    T0->SetElement(1, 3, x0[1]);
    T0->SetElement(2, 3, x0[2]);

    vtkMatrix4x4::Multiply4x4(T0, T2, T3);

    //T3=T1*R1*R2*T0;

    T0->SetElement(0, 3, mps[0][0] - x0[0]);
    T0->SetElement(1, 3, mps[0][1] - x0[1]);
    T0->SetElement(2, 3, mps[0][2] - x0[2]);
    vtkMatrix4x4::Multiply4x4(T0, T3, vtk_T);

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(vtk_T);
    /*
      for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
          cerr << "T(" << i << "," << j << ")=" << vtk_T->GetElement(i,j) << endl;

        }
      }
      //*/
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();

    transformFilter->SetInputData(polydata);
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    // Output EBR mesh for ribSpacingX in project directory
    QString EBRmeshPath = projectDirectory + "/ebr" + QString::number(currentRibSpacing) + ".vtk";
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(EBRmeshPath.toLocal8Bit().data());
    writer->SetInputData(transformFilter->GetOutput());
    writer->Write();

    // Load in EBR mesh
    mesh = mitk::IOUtil::Load<mitk::Surface>(EBRmeshPath.toStdString());
    return mesh;
}

mitk::Surface::Pointer CemrgPower::CalculateAcousticIntensity(mitk::Surface::Pointer endoMesh) {

    /*****************************************************************
    // Calculate acoustic intensity on mesh
    ******************************************************************/
    // Read EBR vtk mesh
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    // Same as output mesh from MapPowerTransmitterToLandmarks
    QString EBRmeshPath = projectDirectory + "/ebr" + QString::number(currentRibSpacing) + ".vtk";

    // ################ NEED if loop to check if file exists! //
    reader->SetFileName(EBRmeshPath.toLocal8Bit().data());
    reader->Update();
    vtkPolyData* ebr_pointSet = reader->GetOutput();

    // MITK load Endo Mesh
    //QString meshPath = projectDirectory + "/transformed-" + QString::number(meshNo) + ".vtk";
    //mitk::Surface::Pointer surf = mitk::IOUtil::Load<mitk::Surface>(meshPath.toStdString());
    //vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();

    vtkSmartPointer<vtkPolyData> endo_polydata = endoMesh->GetVtkPolyData();

    // Read transmitter
    double xaxis[3], yaxis[3], zaxis[3];
    double x0[3], x1[3], y0[3], y1[3], trans_loc[3];

    // Get points from transformed EBR template (Hardcoded)
    //vtkPointSet* ebr_pointSet = transformFilter->GetOutput();
    //xaxis=trans[109129]-trans[70672]
    //yaxis= trans[70672] -trans[62114]

    ebr_pointSet->GetNumberOfPoints();
    //cerr << "Number of points =" << ebr_pointSet->GetNumberOfPoints() << endl;
    //int corner[4]={100247, 72364, 72063, 112175};
    ebr_pointSet->GetPoint(100247, x0);
    ebr_pointSet->GetPoint(72364, x1);
    ebr_pointSet->GetPoint(72063, y0);
    ebr_pointSet->GetPoint(112175, y1);
    //*
    trans_loc[0] = (x0[0] + x1[0] + y0[0] + y1[0]) / 4.0;
    trans_loc[1] = (x0[1] + x1[1] + y0[1] + y1[1]) / 4.0;
    trans_loc[2] = (x0[2] + x1[2] + y0[2] + y1[2]) / 4.0;

    xaxis[0] = x1[0] - x0[0];
    xaxis[1] = x1[1] - x0[1];
    xaxis[2] = x1[2] - x0[2];
    //cerr << "xaxis =" << xaxis(0) << " " << xaxis(1) << " " << xaxis(2) << endl;
    CemrgPower::normalise(xaxis);
    //cerr << "xaxis =" << xaxis(0) << " " << xaxis(1) << " " << xaxis(2) << endl;

    yaxis[0] = x0[0] - y0[0];
    yaxis[1] = x0[1] - y0[1];
    yaxis[2] = x0[2] - y0[2];
    //cerr << "yaxis =" << yaxis(0) << " " << yaxis(1) << " " << yaxis(2) << endl;
    CemrgPower::normalise(yaxis);
    //cerr << "yaxis =" << yaxis(0) << " " << yaxis(1) << " " << yaxis(2) << endl;

    // Find the normal axis to the face of the power transmitter
    //zaxis=yaxis.cross(xaxis);
    //cerr << "y.cross(x) =" << zaxis(0) << " " << zaxis(1) << " " << zaxis(2) << endl;
    CemrgPower::crossProduct(xaxis, yaxis, zaxis);
    //zaxis=xaxis.cross(yaxis);
    //cerr << "x.cross(y) =" << zaxis(0) << " " << zaxis(1) << " " << zaxis(2) << endl;
    CemrgPower::normalise(zaxis);
    //cerr << "zaxis =" << zaxis[0] << " " << zaxis[1] << " " << zaxis[2] << endl;
    // Initialise R;
    vtkSmartPointer<vtkMatrix3x3> R = vtkSmartPointer<vtkMatrix3x3>::New();
    //Eigen::Matrix3f R=Eigen::Matrix3f::Identity();
    vtkSmartPointer<vtkMatrix4x4> R1 = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> vtk_T = vtkSmartPointer<vtkMatrix4x4>::New();

    //Eigen::Matrix4f R1, T;
    // Find rotation matrix to align normal axis of the power transmitter to z-axis
    CemrgPower::fcn_RotationToUnity(zaxis, R);

    R1->SetElement(0, 0, R->GetElement(0, 0));
    R1->SetElement(0, 1, R->GetElement(0, 1));
    R1->SetElement(0, 2, R->GetElement(0, 2));
    R1->SetElement(1, 0, R->GetElement(1, 0));
    R1->SetElement(1, 1, R->GetElement(1, 1));
    R1->SetElement(1, 2, R->GetElement(1, 2));
    R1->SetElement(2, 0, R->GetElement(2, 0));
    R1->SetElement(2, 1, R->GetElement(2, 1));
    R1->SetElement(2, 2, R->GetElement(2, 2));

    // Initialise Transformation matrices;
    vtkSmartPointer<vtkMatrix4x4> T0 = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> T1 = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkMatrix4x4> T2 = vtkSmartPointer<vtkMatrix4x4>::New();
    // Apply translation
    T0->SetElement(0, 3, -trans_loc[0]);
    T0->SetElement(1, 3, -trans_loc[1]);
    T0->SetElement(2, 3, -trans_loc[2]);

    // Scale by 1/1000 mm->m
    T1->SetElement(0, 0, 0.001);
    T1->SetElement(1, 1, 0.001);
    T1->SetElement(2, 2, -0.001);

    // Get transformation to endo mesh
    vtkMatrix4x4::Multiply4x4(R1, T0, T2);
    vtkMatrix4x4::Multiply4x4(T1, T2, vtk_T);

    // Apply transformation to endo mesh so that z-axis is aligned with direction of power beam
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(vtk_T);

    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(endo_polydata);
    transformFilter->SetTransform(transform);
    transformFilter->Update();
    // Get transformed points
    vtkPointSet* endo_pointSet = transformFilter->GetOutput();

    //vtkDataArray *vertex=pointSet->GetPoints()->GetData();

    // Matlab Example Computations of Acoustic Intensity
    // P. Willis EBR Systems Aug 21, 2018
    float f, v, atten, L, A, l;
    f = 921.25E3; // hertz
    v = 1560; // m/sec
    atten = -0.3; //db/(MHz*sec)
    L = 0.9487E-3;
    A = 8 * 24 * pow(L, 2.0);
    l = v / f;

    // Add intensities to each point
    vtkSmartPointer<vtkFloatArray> Intensity = vtkSmartPointer<vtkFloatArray>::New();
    Intensity->SetNumberOfComponents(1);
    Intensity->SetName("Intensity");

    // Get the ids for points where the intensity is > 0
    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);

    std::vector<int> newMeshPoints;
    int size_ids = 0;
    // Calculate the power intensity for each point
    for (int i = 0; i < endo_polydata->GetNumberOfPoints(); i++) {
        double rx[3];
        endo_pointSet->GetPoint(i, rx);

        float d = sqrt(pow(rx[0], 2.0) + pow(rx[1], 2.0) + pow(rx[2], 2.0)); // magnitude of r vector
        float tempx = (M_PI * L * rx[0]) / (l * d);
        float tempy = (M_PI * L * rx[1]) / (l * d);
        float tempI = (A / (pow(l, 2.0) * pow(d, 2.0))) * pow(10.0, (d * atten * f) / 1E5) * sinc(tempx) * sinc(tempy);

        if (tempI >= 0.0) {
            newMeshPoints.push_back(i);
            ids->InsertNextValue(i);
            size_ids++;
            Intensity->InsertNextValue(tempI);
        } else {
            Intensity->InsertNextValue(0.0);
        }
    }
    // Append intensity scalar to input endo mesh
    // endo_polydata->GetPointData()->AddArray(Intensity);
    endo_polydata->GetPointData()->SetScalars(Intensity);
    // sed -i s/'FIELD FieldData 2'/'SCALARS scalars float'/g ~/Dropbox/Work/EBR/data/$case/6_pt/rib/endo${ribspace}.vtk
    // sed -i s/'Intensity [0-9]* [0-9]* double'/'LOOKUP_TABLE default'/g ~/Dropbox/Work/EBR/data/$case/6_pt/rib/endo${ribspace}.vtk

    // So we want to extract cells using points newMeshPoints
    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::POINT);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);
    selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);

    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputConnection(0, reader->GetOutputPort());
    extractSelection->SetInputData(0, endo_polydata);
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    // Write out new endo mesh only where Intensity>0
    QString outEndoMeshPath = projectDirectory + "/endo" + QString::number(currentRibSpacing) + ".vtk";
    vtkSmartPointer<vtkUnstructuredGridWriter> writer2 = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer2->SetFileName(outEndoMeshPath.toLocal8Bit().data());
    writer2->SetInputData(extractSelection->GetOutput());
    writer2->Write();

    // Load in Endo mesh
    mitk::Surface::Pointer mesh;
    mesh = mitk::IOUtil::Load<mitk::Surface>(outEndoMeshPath.toStdString());

    cout << "TESTESTETETET" << mesh->GetVtkPolyData()->GetNumberOfCells() << endl;
    cerr << "output Endo mesh filename in CemrgPower.cpp = " << outEndoMeshPath.toStdString() << endl;
    return mesh;
}

mitk::Surface::Pointer CemrgPower::ReferenceAHA(
    mitk::PointSet::Pointer lmNode, mitk::Surface::Pointer refSurface) {

    //Read the mesh data
    vtkSmartPointer<vtkPolyData> pd = refSurface->GetVtkPolyData();

    //Prepare landmarks
    std::vector<mitk::Point3D> LandMarks;
    if (lmNode.IsNull())
        return refSurface;
    for (mitk::PointSet::PointsIterator it = lmNode->Begin(); it != lmNode->End(); ++it) {
        mitk::Point3D point;
        point.SetElement(0, it.Value().GetElement(0));
        point.SetElement(1, it.Value().GetElement(1));
        point.SetElement(2, it.Value().GetElement(2));
        LandMarks.push_back(point);
    }//for

    // std::vector<mitk::Point3D> LandMarks = ConvertMPS(lmNode);
    // Siemens 4 LM = [apex, basecenter, RV1, RV2]; Siemens 7 LM = [apex, baseMV1, baseMV2, baseMV3, RV1, RV2, apex (just to make it 7 and different from manual LM)]
    // Manual 6 LM = [apex, MV1, MV2, MV3, RV1, RV2]
    if (LandMarks.size() != 4 && LandMarks.size() != 6 && LandMarks.size() != 7) return refSurface;
    mitk::Point3D centre, RIV1, RIV2, APEX, MIV1, MIV2, MIV3;
    if (LandMarks.size() >= 6) {
        APEX = LandMarks.at(0);
        MIV1 = LandMarks.at(1);
        MIV2 = LandMarks.at(2);
        MIV3 = LandMarks.at(3);
        RIV1 = LandMarks.at(4);
        RIV2 = LandMarks.at(5);
        //Calcaulte a circle through the mitral valve points
        centre = Circlefit3d(ZeroPoint(APEX, MIV1), ZeroPoint(APEX, MIV2), ZeroPoint(APEX, MIV3));
        // Zero all points relative to apex
    } else if (LandMarks.size() == 4) {
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

    //Rotate mesh to new frame
    RotateVTKMesh(rotationMat, refSurface);
    //Angle RV cusp 2
    double RVangle1 = atan2(RIV1.GetElement(1), RIV1.GetElement(0));
    double RVangle2 = atan2(RIV2.GetElement(1), RIV2.GetElement(0));
    double appendAngle;
    //Ass