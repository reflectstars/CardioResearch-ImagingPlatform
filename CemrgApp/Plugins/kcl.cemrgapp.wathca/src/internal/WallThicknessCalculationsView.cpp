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
 * Morphological Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>
#include <berryFileEditorInput.h>

// Qmitk
#include <mitkImage.h>
#include <QmitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkBoundingObject.h>
#include <mitkCuboid.h>
#include <mitkAffineImageCropperInteractor.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkUnstructuredGrid.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include "kcl_cemrgapp_wathca_Activator.h"
#include "WallThicknessCalculationsView.h"
#include "WallThicknessCalculationsClipperView.h"

// Micro services
#include <usModuleRegistry.h>

#ifdef _WIN32
// _WIN32 = we're in windows
#include <winsock2.h>
#else
// or linux/mac
#include <arpa/inet.h>
#endif

// VTK
#include <vtkFieldData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>

// ITK
#include <itkAddImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

// CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>

const std::string WallThicknessCalculationsView::VIEW_ID = "org.mitk.views.wathcaview";

WallThicknessCalculationsView::WallThicknessCalculationsView(){
    this->fileName = "";
    this->directory = "";
}

void WallThicknessCalculationsView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(ClipperPV()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(MorphologyAnalysis()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(ThicknessAnalysis()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);
    m_Controls.button_2_4->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    connect(m_Controls.button_2_4, SIGNAL(clicked()), this, SLOT(ApplyFilter()));
    m_Controls.button_3_1->setVisible(false);
    m_Controls.button_3_2->setVisible(false);
    m_Controls.button_3_3->setVisible(false);
    connect(m_Controls.button_3_1, SIGNAL(clicked()), this, SLOT(SelectROI()));
    connect(m_Controls.button_3_2, SIGNAL(clicked()), this, SLOT(SelectLandmarks()));
    connect(m_Controls.button_3_3, SIGNAL(clicked()), this, SLOT(CombineSegs()));
    m_Controls.button_6_1->setVisible(false);
    m_Controls.button_6_2->setVisible(false);
    connect(m_Controls.button_6_1, SIGNAL(clicked()), this, SLOT(ConvertNRRD()));
    connect(m_Controls.button_6_2, SIGNAL(clicked()), this, SLOT(ThicknessCalculator()));
}

void WallThicknessCalculationsView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void WallThicknessCalculationsView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void WallThicknessCalculationsView::Lo