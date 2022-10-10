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

void WallThicknessCalculationsView::LoadDICOM() {

    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void WallThicknessCalculationsView::ProcessIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
        m_Controls.button_2_3->setVisible(false);
        m_Controls.button_2_4->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
        m_Controls.button_2_3->setVisible(true);
        m_Controls.button_2_4->setVisible(true);
    }
}

void WallThicknessCalculationsView::ConvertNII() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() < 1) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please load and select images from the Data Manager before starting this step!");
        return;
    }//_if

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Generic Conversion to nii
    int ctr = 0;
    QString path;
    bool resampleImage = false, reorientToRAI = true;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(nodes.size());
    foreach (mitk::DataNode::Pointer node, nodes) {
        path = directory + "/dcm-" + QString::number(ctr++) + ".nii";
        bool successfulNitfi = CemrgCommonUtils::ConvertToNifti(node->GetData(), path, resampleImage, reorientToRAI);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(node);
        } else {
            mitk::ProgressBar::GetInstance()->Progress(nodes.size());
            return;
        }
        mitk::ProgressBar::GetInstance()->Progress();
    }//_for
    nodes.clear();
    this->BusyCursorOff();

    //Load first item
    ctr = 0;
    path = directory + "/dcm-" + QString::number(ctr) + ".nii";
    mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void WallThicknessCalculationsView::CropIMGS() {

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
                        QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
            if (directory.isEmpty() || directory.simplified().contains(" ")) {
                QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
                directory = QString();
                return;
            }//_if
        }

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

void WallThicknessCalculationsView::ResampIMGS() {

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
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Find the selected node
    QString path;
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            bool ok;
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

            }//_if
        } else
            return;
    } else
        return;
}

void WallThicknessCalculationsView::ApplyFilter() {
}

void WallThicknessCalculationsView::SegmentIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_3_1->isVisible()) {
        m_Controls.button_3_1->setVisible(false);
        m_Controls.button_3_2->setVisible(false);
        m_Controls.button_3_3->setVisible(false);
        return;
    } else {
        m_Controls.button_3_1->setVisible(true);
        m_Controls.button_3_2->setVisible(true);
        m_Controls.button_3_3->setVisible(true);
    }//_if

    int reply = QMessageBox::question(
                NULL, "Question", "Do you have a segmentation to load?",
                QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString path = QFileDialog::getOpenFileName(
                    NULL, "Open Segmentation file", mitk::IOUtil::GetProgramPath().c_str(),
                    QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty())
            return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

        //Restore image name
        QFileInfo fullPathInfo(path);
        fileName = fullPathInfo.fileName();
        directory = fullPathInfo.dir().absolutePath();

    } else {
        //Show the plugin
        QMessageBox::information(
                    NULL, "Attention",
                    "After finalising the segmentation, save the image and load back using this same button!");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

void WallThicknessCalculationsView::SelectROI() {
}

void WallThicknessCalculationsView::SelectLandmarks() {
}

void WallThicknessCalculationsView::CombineSegs() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select both segmentations from the Data Manager to create a surface for clipping!");
        return;
    }//_if

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }//_if

    //Find the selected nodes
    mitk::DataNode::Pointer segNode_0 = nodes.at(0);
    mitk::BaseData::Pointer data_0 = segNode_0->GetData();
    mitk::DataNode::Pointer segNode_1 = nodes.at(1);
    mitk::BaseData::Pointer data_1 = segNode_1->GetData();
    if (data_0 && data_1) {
        //Test if this data item is an image
        mitk::Image::Pointer image_0 = dynamic_cast<mitk::Image*>(data_0.GetPointer());
        mitk::Image::Pointer image_1 = dynamic_cast<mitk::Image*>(data_1.GetPointer());
        if (image_0 && image_1) {

            //Add images
            typedef itk::Image<short, 3> ImageType;
            typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddFilterType;
            QString path = directory + "/segmentation.nii";
            //Cast seg to ITK format
            ImageType::Pointer itkImage_0 = ImageType::New();
            ImageType::Pointer itkImage_1 = ImageType::New();
            CastToItkImage(image_0, itkImage_0);
            CastToItkImage(image_1, itkImage_1);
            AddFilterType::Pointer addFilter = AddFilterType::New();
            addFilter->SetInput1(itkImage_0);
            addFilter->SetInput2(itkImage_1);
            addFilter->UpdateLargestPossibleRegion();

            //Save result
            mitk::Image::Pointer resImage = mitk::ImportItkImage(addFilter->GetOutput())->Clone();
            mitk::IOUtil::Save(resImage, path.toStdString());
            CemrgCommonUtils::AddToStorage(resImage, "segmentation", this->GetDataStorage());
            fileName = "segmentation.nii";

        } else
            return;
    } else
        return;
}

void WallThicknessCalculationsView::ClipperPV() {

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    WallThicknessCalculationsClipperView::SetDirectoryFile(directory, fileName);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.wathcaclipperview");
}

void WallThicknessCalculationsView::MorphologyAnalysis() {

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }//_if

    try {

        QString pathAnalytic = directory + "/AnalyticBloodpool.nii";
        mitk::Image::Pointer analyticImage = mitk::IOUtil::Load<mitk::Image>(pathAnalytic.toStdString());
        QString pathCropped = directory + "/PVeinsCroppedImage.nii";
        mitk::Image::Pointer croppedPVImage = mitk::IOUtil::Load<mitk::Image>(pathCropped.toStdString());

        if (analyticImage && croppedPVImage) {

            //Loop through labelled image
            typedef itk::Image<short, 3> ImageType;
            typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
            ImageType::Pointer analyticItkImage = ImageType::New();
            CastToItkImage(analyticImage, analyticItkImage);
            ItType itLbl(analyticItkImage, analyticItkImage->GetRequestedRegion());
            for (itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itLbl) {
                if ((int)itLbl.Get() == APPENDAGECUT || (int)itLbl.Get() == APPENDAGEUNCUT) {
                    itLbl.Set(0);
                }//_if
            }//_for

            //Relabel the components to separate bloodpool and appendage
            typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
            ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
            connected->SetInput(analyticItkImage);
            connected->Update();
            typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
            RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
            relabeler->SetInput(connected->GetOutput());
            relabeler->Update();

            //Keep the selected labels
            typedef itk::LabelObject<short, 3> LabelObjectType;
            typedef itk::LabelMap<LabelObjectType> LabelMapType;
            typedef itk::LabelImageToLabelMapFilter< ImageType, LabelMapType > LabelImageToLabelMapFilterType;
            LabelImageToLabelMapFilterType::Pointer labelMapConverter = LabelImageToLabelMapFilterType::New();
            labelMapConverter->SetInput(relabeler->GetOutput());
            labelMapConverter->SetBackgroundValue(0);
            typedef itk::LabelSelectionLabelMapFilter<LabelMapType> SelectorType;
            SelectorType::Pointer selector = SelectorType::New();
            selector->SetInput(labelMapConverter->GetOutput());
            selector->SetLabel(2);

            //Import to MITK image
            typedef itk::LabelMapToLabelImageFilter<LabelMapType, ImageType> LabelMapToLabelImageFilterType;
            LabelMapToLabelImageFilterType::Pointer labelImageConverter = LabelMapToLabelImageFilterType::New();
            labelImageConverter->SetInput(selector->GetOutput(0));
            labelImageConverter->Update();
            mitk::Image::Pointer ap = mitk::ImportItkImage(labelImageConverter->GetOutput());
            mitk::Image::Pointer bp = mitk::IOUtil::Load<mitk::Image>(directory.toStdString() + "/PVeinsCroppedImage.nii");

            //Ask for user input to set the parameters
            QDialog* inputs = new QDialog(0,0);
            m_UIMeshing.setupUi(inputs);
            connect(m_UIMeshing.buttonBox, &QDialogButtonBox::accepted, inputs, &QDialog::accept);
            connect(m_UIMeshing.buttonBox, &QDialogButtonBox::rejected, inputs, &QDialog::reject);
            int dialogCode = inputs->exec();

            //Act on dialog return code
            if (dialogCode == QDialog::Accepted) {

                bool ok1, ok2, ok3, ok4;
                float th = m_UIMeshing.lineEdit_1->text().toFloat(&ok1);
                float bl = m_UIMeshing.lineEdit_2->text().toFloat(&ok2);
                int smth = m_UIMeshing.lineEdit_3->text().toInt(&ok3);
                float ds = m_UIMeshing.lineEdit_4->text().toFloat(&ok4);

                //Set default values
                if (!ok1 || !ok2 || !ok3 || !ok4)
             