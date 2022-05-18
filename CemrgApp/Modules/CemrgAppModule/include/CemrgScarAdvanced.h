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
 * Scar Advanced Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * jose.solislemus@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgScarAdvanced_h
#define CemrgScarAdvanced_h

// Qmitk
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <vtkFloatArray.h>
#include <MitkCemrgAppModuleExports.h>

// VTK
#include <vtkAppendFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkConnectivityFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkCellLocator.h>
#include "vtkPolyData.h"
#include <vtkCellArray.h>
#include "vtkGenericCell.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include <vtkCellData.h>
#include <vtkPointLocator.h>
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataReader.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkThreshold.h>
#include <vtkThresholdPoints.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkTriangle.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkSelectPolyData.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkStructuredPoints.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLineSource.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointPicker.h>

// C++ Standard
#include <string>
#include <sstream>

class MITKCEMRGAPPMODULE_EXPORT CemrgScarAdvanced {

public:

    bool _debugScarAdvanced;

    int _neighbourhood_size;
    int _run_count;
    double _fill_threshold;
    double _max_scalar;
    bool _weightedcorridor;
    std::string _fileOutName;
    std::string _outPath;
    std::string _prefix;
    std::string _leftrightpre;

    double fi1_largestSurfaceArea, fi1_scarScore;
    double fi2_percentage, fi2_largestSurfaceArea, fi2_corridorSurfaceArea;
    int fi2_connectedAreasTotal;
    double fi3_preScarScoreSimple, fi3_postScarScoreSimple;
    double fi3_totalPoints, fi3_emptyPoints, fi3_healthy, fi3_preScar, fi3_postScar, fi3_overlapScar;
    std::string fi1_fname, fi2_fname, fi3_fname;

    std::vector<std::pair<int, int> > _visited_point_list; // stores the neighbours around a point
    std::vector<vtkSmartPointer<vtkPolyData> > _paths; // container to store shortest paths between points
    std::vector<vtkSmartPointer<vtkPolyDataMapper> > _pathMappers;
    std::vector<vtkSmartPointer<vtkActor> > _actors;

    vtkSmartPointer<vtkCellPicker> _cell_picker;
    vtkSmartPointer<vtkPointPicker> _point_picker;
    vtkSmartPointer<vtkPolyData> _SourcePolyData;

    vtkSmartPointer<vtkPolyData> _source;
    vtkSmartPointer<vtkPolyData> _target;

    std::vector<vtkSmartPointer<vtkDijkstraGraphGeodesicPath> > _shortestPaths;
    std::vector<int> _pointidarray;
    std::vector<int> _corridoridarray;

    // Getters and setters
    inline bool IsDebug() { return _debugScarAdvanced; };

    inline bool IsWeighted