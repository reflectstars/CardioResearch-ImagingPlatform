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
 * Anatomical Measurements
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
#include <mitkCuboid.h>
#include <mitkDataNode.h>
#include <mitkProgressBar.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkAffineImageCropperInteractor.h>
#include "kcl_cemrgapp_mmeasurement_Activator.h"
#include "MmeasurementView.h"

// Micro services
#include <usModuleRegistry.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QSignalMapper>
#include <QInputDialog>

// C++ Standard
#include <numeric>

// CemrgAppModule
#include <CemrgMeasure.h>
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>

const std::string MmeasurementView::VIEW_ID = "org.mitk.views.motionmeasurement";

MmeasurementView::MmeasurementView() {
    this->timePoints = 0;
    this->smoothness = 0;
    this->directory = "";
}

void MmeasurementView::CreateQtPartControl(QWidget *parent) {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(TrackingButton()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(SelectLandmark()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(WriteFileButton()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(CalcDistButton()));
    connect(m_Controls.button_7, SIGNAL(clicked()), this, SLOT(CalcPeriButton()));
    connect(m_Controls.button_8, SIGNAL(clicked()), this, SLOT(CalcAreaButton()));
    connect(m_Controls.button_9, SIGNAL(clicked()), this, SLOT(FindCentre()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Processing buttons
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropinIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);

    //Tracking buttons
    connect(m_Controls.button_3_1, SIGNAL(clicked()), this, SLOT(Tracking()));
    connect(m_Controls.button_3_2, SIGNAL(clicked()), this, SLOT(Applying()));
    //Set visibility of buttons
    m_Controls.button_3_1->setVisible(false);
    m_Controls.button_3_2->setVisible(false);

    //Temporal resolution
    timePoints = 0;
    smoothness = 0;
}

void MmeasurementView::SetFocus() {
    m_Controls.button_1->setFocus();
}

void MmeasurementView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void MmeasurementView::LoadDICOM() {
    // Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void MmeasurementView::ProcessIMGS() {
    // Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
        m_Contro