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
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void EASIView::ProcessIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
        m_Controls.button_2_3->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
        m_Controls.button_2_3->setVisible(true);
    }
}

void EASIView::ConvertNII() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 10) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please load and select 10 images from the Data Manager before starting this step!");
        return;
    }//_if

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
            NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Order dicoms based on their cycle stages
    std::vector<int> indexNodes;
    std::string seriesDescription;
    foreach (mitk::DataNode::Pointer node, nodes) {
        node->GetData()->GetPropertyList()->GetStringProperty("dicom.series.SeriesDescription", seriesDescription);
        if (seriesDescription.find("90.0%") != seriesDescription.npos) indexNodes.push_back(9);
        else if (seriesDescription.find("80.0%") != seriesDescription.npos) indexNodes.push_back(8);
        else if (seriesDescription.find("70.0%") != seriesDescription.npos) indexNodes.push_back(7);
        else if (seriesDescription.find("60.0%") != seriesDescription.npos) indexNodes.push_back(6);
        else if (seriesDescription.find("50.0%") != seriesDescription.npos) indexNodes.push_back(5);
        else if (seriesDescription.find("40.0%") != seriesDescription.npos) indexNodes.push_back(4);
        else if (seriesDescription.find("30.0%") != seriesDescription.npos) indexNodes.push_back(3);
        else if (seriesDescription.find("20.0%") != seriesDescription.npos) indexNodes.push_back(2);
        else if (seriesDescription.find("10.0%") != seriesDescription.npos) indexNodes.push_back(1);
        else if (seriesDescription.find("0.0%") != seriesDescription.npos) indexNodes.push_back(0);
    }//_for

    //Sort indexes based on comparing values
    std::vector<int> index(indexNodes.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](int i1, int i2) {return indexNodes[i1] < indexNodes[i2]; });
    //Warning for cases when order is not found
    size_t length1 = nodes.size();
    size_t length2 = 