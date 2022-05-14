
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
 * Simple Common Utilities
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgCommonUtils_h
#define CemrgCommonUtils_h

#include <MitkCemrgAppModuleExports.h>
#include <mitkImage.h>
#include <mitkBoundingObject.h>
#include <mitkDataNode.h>
#include <mitkDataStorage.h>
#include <QString>

class MITKCEMRGAPPMODULE_EXPORT CemrgCommonUtils {

public:

    //Cropping Utils
    static mitk::Image::Pointer CropImage();
    static void SetImageToCut(mitk::Image::Pointer imageToCut);
    static void SetCuttingCube(mitk::BoundingObject::Pointer cuttingCube);
    static void SetImageNode(mitk::DataNode::Pointer imageNode);
    static void SetCuttingNode(mitk::DataNode::Pointer cuttingNode);
    static mitk::DataNode::Pointer GetImageNode();
    static mitk::DataNode::Pointer GetCuttingNode();

    //Sampling Utils
    static mitk::Image::Pointer Downsample(mitk::Image::Pointer image, int factor);
    static mitk::Image::Pointer IsoImageResampleReorient(mitk::Image::Pointer image, bool resample = false, bool reorientToRAI = false);
    static mitk::Image::Pointer IsoImageResampleReorient(QString imPath, bool resample = false, bool reorientToRAI = false);

    // Image Analysis Utils
    static void SetSegmentationEdgesToZero(mitk::Image::Pointer image, QString outPath = "");

    //Nifti Conversion Utils
    static bool ConvertToNifti(mitk::BaseData::Pointer oneNode, QString path2file, bool resample = false, bool reorient = false);
    static void RoundPixelValues(QString pathToImage, QString outputPath = "");

    // static void RoundPointDataValues(vtkSmartPointer<vtkPolyData> pd);

    //Mesh Utils
    static mitk::Surface::Pointer LoadVTKMesh(std::string path);
    static mitk::Surface::Pointer ExtractSurfaceFromSegmentation(mitk::Image::Pointer image, double thresh = 0.5, double blur = 0.8, double smoothIterations = 3, double decimation = 0.5);
    static mitk::Surface::Pointer ClipWithSphere(mitk::Surface::Pointer surface, double x_c, double y_c, double z_c, double radius, QString saveToPath = "");
    static void FlipXYPlane(mitk::Surface::Pointer surf, QString dir, QString vtkname = "segmentation.vtk");
    static QString M3dlibParamFileGenerator(QString dir, QString filename = "param-template.par", QString thicknessCalc = "0");
    static bool ConvertToCarto(std::string vtkPath, std::vector<double>, double, double, int, bool);
    static void CalculatePolyDataNormals(vtkSmartPointer<vtkPolyData>& pd, bool celldata = true);
    static void FillHoles(mitk::Surface::Pointer surf, QString dir = "", QString vtkname = "");

    //Tracking Utils
    static void MotionTrackingReport(QString directory, int timePoints);

    //Generic
    static mitk::DataNode::Pointer AddToStorage(mitk::BaseData* data, std::string nodeName, mitk::DataStorage::Pointer ds, bool init = true);

    //Carp Utils
    static void OriginalCoordinates(QString imagePath, QString pointPath, QString outputPath, double scaling = 1000);
    static void CalculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath);
    static void RegionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath);
    static void NormaliseFibreFiles(QString fibresPath, QString outputPath);
    static void RectifyFileValues(QString pathToFile, double minVal = 0.0, double maxVal = 1.0);
    static int GetTotalFromCarpFile(QString pathToFile, bool totalAtTop = true);
    static std::vector<double> ReadScalarField(QString pathToFile);
    static void CarpToVtk(QString elemPath, QString ptsPath, QString outputPath, bool saveRegionlabels = true);
    static void AppendScalarFieldToVtk(QString vtkPath, QString fieldName, QString typeData, std::vector<double> field, bool setHeader = true);
    static void AppendVectorFieldToVtk(QString vtkPath, QString fieldName, QString typeData, std::vector<double> field, bool setHeader = true);

private:

    //Cropping Utils
    static mitk::Image::Pointer imageToCut;
    static mitk::BoundingObject::Pointer cuttingCube;
    static mitk::DataNode::Pointer imageNode;
    static mitk::DataNode::Pointer cuttingNode;
};

#endif // CemrgCommonUtils_h