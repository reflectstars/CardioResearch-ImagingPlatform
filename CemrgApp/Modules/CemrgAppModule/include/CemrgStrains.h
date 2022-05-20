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

#ifndef CemrgStrains_h
#define CemrgStrains_h

#include <mitkSurface.h>
#include <vtkCell.h>
#include <vtkFloatArray.h>
#include <MitkCemrgAppModuleExports.h>

class MITKCEMRGAPPMODULE_EXPORT CemrgStrains {

public:

    CemrgStrains();
    CemrgStrains(QString dir, int refMeshNo);
    ~CemrgStrains();

    double CalculateGlobalSqzPlot(int meshNo);
    std::vector<double> CalculateSqzPlot(int meshNo);
    std::vector<double> CalculateStrainsPlot(int meshNo, mitk::DataNode::Pointer lmNode, int flag);
    double CalculateSDI(std::vector<std::vector<double>> valueVectors, int cycleLengths, int noFrames);

    std::vector<mitk::Surface::Pointer> ReferenceGuideLines(mitk::DataNode::Pointer lmNode);
    mitk::Surface::Pointer ReferenceAHA(mitk::DataNode::Pointer lmNode, int segRatios[], bool pacingSite);
    mitk::Surface::Pointer FlattenedAHA();
    vtkSmartPointer<vtkFloatArray> GetFlatSurfScalars() const;
    std::vector<float> GetAHAColour(int label);

protected:

    mitk::Point3D ZeroPoint(mitk::Point3D apex, mitk::Point3D point);
    void ZeroVTKMesh(mitk::Point3D apex, mitk::Surface: