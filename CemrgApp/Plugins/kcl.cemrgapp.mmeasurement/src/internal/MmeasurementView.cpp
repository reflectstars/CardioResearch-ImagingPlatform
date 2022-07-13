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
        m_Controls.button_2_3->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
        m_Controls.button_2_3->setVisible(true);
    }
}

void MmeasurementView::ConvertNII() {
    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != timePoints) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please load and select all images from the Data Manager before starting this step!");
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
    size_t length2 = indexNodes.size();
    if (length1 != length2) {
        QMessageBox::warning(
            NULL, "Attention",
            "Cannot find the order of images automatically. Revert to user order and selections in the data manager!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString path;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        path = directory + "/dcm-" + QString::number(ctr++) + ".nii";
        bool successfulNitfi = CemrgCommonUtils::ConvertToNifti(nodes.at(idx)->GetData(), path);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(nodes.at(idx));
        } else {
            mitk::ProgressBar::GetInstance()->Progress(index.size());
            return;
        }
        mitk::ProgressBar::GetInstance()->Progress();
    }//for
    nodes.clear();
    this->BusyCursorOff();

    //Load first item
    ctr = 0;
    path = directory + "/dcm-" + QString::number(ctr) + ".nii";
    mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void MmeasurementView::CropinIMGS() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please select an image from the Data Manager to perform cropping!");
        return;
    }//_if

    //Check to cut now or not
    if (m_Controls.button_2_2->text() == QString::fromStdString("Are you done?")) {

        QString path;
        //Ask the user for a dir to locate data
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

        //Check for temporal resolution
        bool ok = true;
        if (timePoints == 0)
            timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
        if (!ok) {
            QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
            timePoints = 0;
            return;
        }//_if

        //Cut selected image
        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
        mitk::Image::Pointer outputImage = CemrgCommonUtils::CropImage();
        path = directory + "/" + CemrgCommonUtils::GetImageNode()->GetName().c_str() + ".nii";
        mitk::IOUtil::Save(outputImage, path.toStdString());
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();

        //Update datastorage
        CemrgCommonUtils::AddToStorage(outputImage, CemrgCommonUtils::GetImageNode()->GetName(), this->GetDataStorage());
        this->GetDataStorage()->Remove(CemrgCommonUtils::GetImageNode());
        this->GetDataStorage()->Remove(CemrgCommonUtils::GetCuttingNode());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

        //Cut rest of images
        int reply = QMessageBox::question(
            NULL, "Question", "Would you like to automate cropping of other images in the cycle?",
            QMessageBox::Yes, QMessageBox::No);
        if (reply == QMessageBox::Yes) {

            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(timePoints - 1);

            for (int i = 1; i < timePoints; i++) {

                mitk::Image::Pointer inputImage;
                path = directory + "/dcm-" + QString::number(i) + ".nii";
                try {
                    inputImage = dynamic_cast<mitk::Image*>(mitk::IOUtil::Load(path.toStdString()).front().GetPointer());
                } catch (const std::exception&) {
                    mitk::ProgressBar::GetInstance()->Progress();
                    continue;
                }//_try

                //Setup cropper
                CemrgCommonUtils::SetImageToCut(inputImage);
                outputImage = CemrgCommonUtils::CropImage();
                mitk::IOUtil::Save(outputImage, path.toStdString());
                mitk::ProgressBar::GetInstance()->Progress();

            }//_for
            this->BusyCursorOff();
        }//_if

        m_Controls.button_2_2->setText("Crop Images");
        return;
    }//_if

    //Prepare cutting cuboid
    mitk::Cuboid::Pointer cuttingObject = mitk::Cuboid::New();
    mitk::DataNode::Pointer cuttingNode = mitk::DataNode::New();
    cuttingNode->SetData(cuttingObject);
    cuttingNode->SetProperty("opacity", mitk::FloatProperty::New(0.4));
    cuttingNode->SetProperty("color", mitk::ColorProperty::New(1.0, 1.0, 0.0));
    cuttingNode->SetProperty("name", mitk::StringProperty::New("Cropper"));
    this->GetDataStorage()->Add(cuttingNode);

    //Mouse interactions
    mitk::AffineImageCropperInteractor::Pointer affineDataInteractor = mitk::AffineImageCropperInteractor::New();
    affineDataInteractor->LoadStateMachine("ClippingPlaneInteraction3D.xml", us::ModuleRegistry::GetModule("MitkDataTypesExt"));
    affineDataInteractor->SetEventConfig("CropperDeformationConfig.xml", us::ModuleRegistry::GetModule("MitkDataTypesExt"));
    affineDataInteractor->SetDataNode(cuttingNode);
    cuttingNode->SetBoolProperty("pickable", true);

    //Fit the cuboid to the image
    mitk::Image::Pointer imageToCut;
    mitk::BoundingObject::Pointer cuttingCube;
    mitk::DataNode::Pointer imageNode = nodes.at(0);
    mitk::BaseData::Pointer data = imageNode->GetData();
    cuttingCube = dynamic_cast<mitk::BoundingObject*>(cuttingNode->GetData());
    if (data) {
        //Test if this data item is an image
        imageToCut = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (imageToCut)
            cuttingCube->FitGeometry(imageToCut->GetGeometry());
        else return;
    } else return;

    //To be used for actual cutting
    CemrgCommonUtils::SetImageToCut(imageToCut);
    CemrgCommonUtils::SetCuttingCube(cuttingCube);
    CemrgCommonUtils::SetImageNode(imageNode);
    CemrgCommonUtils::SetCuttingNode(cuttingNode);
    m_Controls.button_2_2->setText("Are you done?");
}

void MmeasurementView::ResampIMGS() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please select an image from the Data Manager to perform downsampling!");
        return;
    }

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

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if

    //Find the selected node
    QString path;
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {
            int factor = QInputDialog::getInt(NULL, tr("Downsampling"), tr("By factor of:"), 3, 1, 5, 1, &ok);
            if (ok) {
                //Downsample selected image
                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
                mitk::Image::Pointer outputImage = CemrgCommonUtils::Downsample(image, factor);
                path = directory + "/" + imgNode->GetName().c_str() + ".nii";
                mitk::IOUtil::Save(outputImage, path.toStdString());
                mitk::ProgressBar::GetInstance()->Progress();
                this->BusyCursorOff();

                //Update datastorage
                CemrgCommonUtils::AddToStorage(outputImage, imgNode->GetName(), this->GetDataStorage());
                this->GetDataStorage()->Remove(imgNode);
                mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

                //Downsample rest of images
                int reply = QMessageBox::question(
                    NULL, "Question", "Would you like to automate downsampling of other images in the cycle?",
                    QMessageBox::Yes, QMessageBox::No);

                if (reply == QMessageBox::Yes) {

                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(timePoints - 1);
                    for (int i = 1; i < timePoints; i++) {

                        mitk::Image::Pointer inputImage;
                        path = directory + "/dcm-" + QString::number(i) + ".nii";
                        try {
                            inputImage = dynamic_cast<mitk::Image*>(mitk::IOUtil::Load(path.toStdString()).front().GetPointer());
                        } catch (const std::exception&) {
                            mitk::ProgressBar::GetInstance()->Progress();
                            continue;
                        }//_try

                        //Setup sampler
                        outputImage = CemrgCommonUtils::Downsample(inputImage, factor);
                        mitk::IOUtil::Save(outputImage, path.toStdString());
                        mitk::ProgressBar::GetInstance()->Progress();

                    }//_for
                    this->BusyCursorOff();
                }//_if
            }//_if
        } else
            return;
    } else
        return;
}

void MmeasurementView::TrackingButton() {

    //Toggle visibility of buttons
    if (m_Controls.button_3_1->isVisible()) {
        m_Controls.button_3_1->setVisible(false);
        m_Controls.button_3_2->setVisible(false);
    } else {
        m_Controls.button_3_1->setVisible(true);
        m_Controls.button_3_2->setVisible(true);
    }
}

void MmeasurementView::SelectLandmark() {

    //Show the plugin
    QMessageBox::information(NULL, "Attention", "Please select your points in order!");
    this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
}

void MmeasurementView::BrowseT(const QString& buttDir) {

    QString time, para = "";
    QString buttID = buttDir.left(1);
    QString direct = buttDir.right(buttDir.size() - 1);

    //Load target, time and parameter files
    switch (buttID.toInt()) {
    case 1:
        time = QFileDialog::getOpenFileName(
            NULL, "Open text file containing time points of source images",
            direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UITracking.lineEdit_1->setText(time);
        break;
    case 2:
        para = QFileDialog::getOpenFileName(
            NULL, "Open text file containing parameters",
            direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UITracking.lineEdit_2->setText(para);
        break;
    }//_switch
}

void MmeasurementView::BrowseA(const QString& buttDir) {

    QString buttID = buttDir.left(1);
    QString direct = buttDir.right(buttDir.size() - 1);

    //Load input mesh, dofin file
    switch (buttID.toInt()) {
    case 1:
        // QString input = QFileDialog::getOpenFileName(
        //     NULL, "Open the input mesh",
        //     direct, QmitkIOUtil::GetFileOpenFilterString());
        //m_UIApplying.lineEdit_1->setText(input);
        break;
    case 2:
        QString dofin = QFileDialog::getOpenFileName(
            NULL, "Open the transformation file",
            direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UIApplying.lineEdit_3->setText(dofin);
        break;
    }//_switch
}

void MmeasurementView::Tracking() {

    //Ask the user for project directory
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
            NULL, "Open Project Directory to Locate Images", mitk::IOUtil::GetProgramPath().c_str(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0, 0);
    QSignalMapper* signalMapper = new QSignalMapper(this);

    m_UITracking.setupUi(inputs);
    connect(m_UITracking.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UITracking.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_UITracking.pushButton_1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    connect(m_UITracking.pushButton_2, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_UITracking.pushButton_1, "1" + directory);
    signalMapper->setMapping(m_UITracking.pushButton_2, "2" + directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(BrowseT(const QString&)));

    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        QString time = m_UITracking.lineEdit_1->text();
        QString para = m_UITracking.lineEdit_2->text();

        //Checking input files
        if (time.isEmpty() || para.isEmpty())
            QMessageBox::warning(NULL, "Attention", "Reverting to default time or parameter file!");

        if (time.isEmpty()) {
            ofstream file;
            //Absolute path
            QString aPath = QCoreApplication::applicationDirPath() + "/MLib";
            file.open(aPath.toStdString() + "/imgTimes.lst");
            file << directory << "/dcm- .nii" << endl;
            for (int i = 0; i < timePoints; i++)
                file << i << " " << i * 10 << endl;
            file.close();
            time = aPath + "/imgTimes.lst";
        }//_if

        //Commandline execution
        this->BusyCursorOn();
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->ExecuteTracking(directory, time, para);
        QMessageBox::information(NULL, "Attention"