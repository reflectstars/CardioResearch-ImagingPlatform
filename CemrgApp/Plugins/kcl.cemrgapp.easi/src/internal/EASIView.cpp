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
 * Eikonal Activation Simulation (EASI) Plugin for MITK
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
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkBoundingObject.h>
#include <mitkCuboid.h>
#include <mitkEllipsoid.h>
#include <mitkAffineImageCropperInteractor.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkUnstructuredGrid.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImageCast.h>
#include <mitkScaleOperation.h>
#include <mitkInteractionConst.h>
#include <mitkImage.h>
#include "kcl_cemrgapp_easi_Activator.h"
#include "EASIView.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

// VTK
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// CemrgAppModule
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>
// #include "CemrgTests.cpp" // Might have to be accessed directly to a path in the developer's specific workstation

// C++ Standard
#include <usModuleRegistry.h>
#include <numeric>
#include <iostream>

const std::string EASIView::VIEW_ID = "org.mitk.views.easi";

void EASIView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void EASIView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(CreateMesh()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(ActivationSites()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(Simulation()));
    connect(m_Controls.button_l, SIGNAL(clicked()), this, SLOT(LoadMesh()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropinIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    m_Controls.button_5_1->setVisible(false);
    connect(m_Controls.button_5_1, SIGNAL(clicked()), this, SLOT(ConfrmSITE()));
}

void EASIView::OnSelectionChanged(
    berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void EASIView::LoadDICOM() {

    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.di