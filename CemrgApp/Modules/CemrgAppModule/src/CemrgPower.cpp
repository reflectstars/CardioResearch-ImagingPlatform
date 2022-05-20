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

    // Check that there are more than 2