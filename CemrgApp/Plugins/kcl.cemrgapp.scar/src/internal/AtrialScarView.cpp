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
 * Scar Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryFileEditorInput.h>
#include <berryIWorkbenchPage.h>
#include <berryISelectionProvider.h>
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include <mitkImage.h>
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkLookupTableProperty.h>
#include <mitkVtkScalarModeProperty.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include "kcl_cemrgapp_scar_Activator.h"
#include "AtrialScarView.h"
#include "AtrialScarClipperView.h"
#include "ScarCalculationsView.h"

// VTK
#include <vtkPolyData.h>
#include <vtkLookupTable.h>
#include <vtkSphereSource.h>
#include <vtkImageResize.h>
#include <vtkImageChangeInformation.h>
#include <vtkCellDataToPointData.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkDecimatePro.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkTimerLog.h>
#include <vtkPPolyDataNormals.h>

// ITK
#include <itkResampleImageFilter.h>
#include <itkBinaryCrossStructuringElement.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkImageFileWriter.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QStringList>

// C++ Standard
#include <numeric>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>
#include <CemrgCommonUtils.h>

const std::string AtrialScarView::VIEW_ID = "org.mitk.views.scar";

AtrialScarView::AtrialScarView(){
    this->fileName = "";
    this->directory = "";
    this->debugSCARname = "";
    this->alternativeNiftiFolder = "";
}

void AtrialScarView::CreateQtPartControl(QWidget *parent) {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(AnalysisChoice()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(CreateSurf()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(ScarMap()));
    connect(m_Controls.button_7, SIGNAL(clicked()), this, SLOT(Threshold()));
    connect(m_Controls.button_x, SIGNAL(clicked()), this, SLOT(Registration()));
    connect(m_Controls.button_y, SIGNAL(clicked()), this, SLOT(ClipPVeins()));
    connect(m_Controls.button_z, SIGNAL(clicked()), this, SLOT(ClipperMV()));
    connect(m_Controls.button_s, SIGNAL(clicked()), this, SLOT(Sphericity()));
    connect(m_Controls.button_c, SIGNAL(clicked()), this, SLOT(ExtraCalcs()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(ResetMain()));

    //Sub-buttons signals
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_x_1, SIGNAL(clicked()), this, SLOT(Register()));
    connect(m_Controls.button_x_2, SIGNAL(clicked()), this, SLOT(Transform()));
    connect(m_Controls.button_z_1, SIGNAL(clicked()), this, SLOT(SelectLandmarks()));
    connect(m_Controls.button_z_2, SIGNAL(clicked()), this, SLOT(ClipMitralValve()));
    connect(m_Controls.button_deb, SIGNAL(clicked()), this, SLOT(ScarDebug()));

    //Set visibility of buttons
    m_Controls.button_4->setVisible(false);
    m_Controls.button_5->setVisible(false);
    m_Controls.button_6->setVisible(false);
    m_Controls.button_7->setVisible(false);
    m_Controls.button_x->setVisible(false);
    m_Controls.button_y->setVisible(false);
    m_Controls.button_z->setVisible(false);
    m_Controls.button_s->setVisible(false);

    //Set visibility of sub-buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_x_1->setVisible(false);
    m_Controls.button_x_2->setVisible(false);
    m_Controls.button_z_1->setVisible(false);
    m_Controls.button_z_2->setVisible(false);
    m_Controls.button_deb->setVisible(false);
}

void AtrialScarView::SetFocus() {
    m_Controls.button_1->setFocus();
}

void AtrialScarView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void AtrialScarView::LoadDICOM() {

    MITK_INFO << "Ask user about alternative DICOM reader";
    int reply = QMessageBox::question(NULL, "Question", "Use alternative DICOM reader?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString dicomFolder = QFileDialog::getExistingDirectory(NULL, "Open folder with DICOMs.", mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        QString tmpNiftiFolder = cmd->DockerDicom2Nifti(dicomFolder);

        if (tmpNiftiFolder.compare("ERROR_IN_PROCESSING") != 0) {

            // add results in NIIs folder to Data Manager
            MITK_INFO << ("Conversion succesful. Intermediate NII folder: " + tmpNiftiFolder).toStdString();
            QMessageBox::information(NULL, "Information", "Conversion successful, press the Process Images button to continue.");
            QDir niftiFolder(tmpNiftiFolder);
            QStringList niftiFiles = niftiFolder.entryList();

            if (niftiFiles.size() > 0) {

                QString thisFile, path;
                for (int ix = 0; ix < niftiFiles.size(); ix++) {

                    // load here files
                    thisFile = niftiFiles.at(ix);
                    if (thisFile.contains(".nii", Qt::CaseSensitive)) {
                        if (thisFile.contains("lge", Qt::CaseInsensitive) || thisFile.contains("mra", Qt::CaseInsensitive)) {

                            path = niftiFolder.absolutePath() + "/" + thisFile;
                            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
                            std::string key = "dicom.series.SeriesDescription";
                            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), thisFile.left(thisFile.length() - 4).toStdString().c_str());

                        }//_if
                    }//_if
                }//_for

            } else {

                MITK_WARN << "Problem with conversion.";
                QMessageBox::warning(NULL, "Attention", "Problem with alternative conversion. Try MITK Dicom editor?");
                return;

            }//_if
        }//_if

    } else {

        MITK_INFO << "Using MITK DICOM editor";
        QString editor_id = "org.mitk.editors.dicomeditor";
        berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
        this->GetSite()->GetPage()->OpenEditor(input, editor_id);

    }//_if
}

void AtrialScarView::ProcessIMGS() {
    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
    }
}

void AtrialScarView::ConvertNII() {
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        QMessageBox::warning(NULL, "Attention", "Please load and select both LGE and CEMRA images from the Data Manager to convert!");
        return;
    }//_if

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Order dicoms based on their type
    std::vector<int> indexNodes;
    std::vector<std::string> seriesDscrps;
    foreach (mitk::DataNode::Pointer node, nodes) {

        std::string seriesDescription;
        node->GetData()->GetPropertyList()->GetStringProperty("dicom.series.SeriesDescription", seriesDescription);

        if (seriesDescription.find("LGE") != seriesDescription.npos) indexNodes.push_back(0);
        else if (seriesDescription.find("MRA") != seriesDescription.npos) indexNodes.push_back(1);

        //Trim whitespaces
        seriesDescription = QString::fromStdString(seriesDescription).replace(")", "").toStdString();
        seriesDescription = QString::fromStdString(seriesDescription).replace("(", "").toStdString();
        seriesDescription = QString::fromStdString(seriesDescription).simplified().replace(" ", "").toStdString();
        seriesDscrps.push_back(seriesDescription);
    }//_for

    //Sort indexes based on comparing values
    std::vector<int> index(indexNodes.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](int i1, int i2) {return indexNodes[i1] < indexNodes[i2]; });

    //Warning for cases when type is not found
    size_t length1 = nodes.size();
    size_t length2 = indexNodes.size();
    bool test = std::adjacent_find(indexNodes.begin(), indexNodes.end(), std::not_equal_to<int>()) == indexNodes.end();
    if (length1 != length2 || test) {
        QMessageBox::warning(NULL, "Attention",
            "Cannot find the type of images automatically. Revert to user order and selections in the data manager: LGE at the top, then CEMRA at the bottom!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString path, type;
    bool resampleImage = true, reorientToRAI = true;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        type = (ctr == 0) ? "LGE" : "MRA";
        path = directory + "/dcm-" + type + "-" + seriesDscrps.at(idx).c_str() + ".nii";
        bool successfulNitfi = CemrgCommonUtils::ConvertToNifti(nodes.at(idx)->GetData(), path, resampleImage, reorientToRAI);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(nodes.at(idx));
            std::string key = "dicom.series.SeriesDescription";
            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), seriesDscrps.at(idx).c_str());
            ctr++;
        } else {
            mitk::ProgressBar::GetInstance()->Progress(index.size());
            return;
        }//_if
        mitk::ProgressBar::GetInstance()->Progress();
    }//for
    nodes.clear();
    this->BusyCursorOff();

    MITK_INFO << "Loading all items";
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void AtrialScarView::AnalysisChoice() {

    int reply = QMessageBox::question(
        NULL, "Question", "Do you want an automatic analysis?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        MITK_INFO << "Setting up automatic analysis.";
        AutomaticAnalysis();

    } else {

        MITK_INFO << "Setting up manual analysis.";

        //Set visibility of buttons
        m_Controls.button_4->setVisible(true);
        m_Controls.button_5->setVisible(true);
        m_Controls.button_6->setVisible(true);
        m_Controls.button_7->setVisible(true);
        m_Controls.button_x->setVisible(true);
        m_Controls.button_y->setVisible(true);
        m_Controls.button_z->setVisible(true);
        m_Controls.button_s->setVisible(true);
    }//_if
}

void AtrialScarView::AutomaticAnalysis() {

    MITK_INFO << "Performing automatic analysis.";
    MITK_INFO << "============= Automatic segmentation module ====================";

    QString direct, mraPath, lgePath, cnnPath;
    bool debugging = true;

    if (directory.isEmpty()) {
        direct = QFileDialog::getExistingDirectory(
            NULL, "Open Project Directory",
            mitk::IOUtil::GetProgramPath().c_str(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
        directory = direct;
    } else {
        direct = directory;
    }//_dir

    QDirIterator searchit(direct, QDirIterator::Subdirectories);

    MITK_INFO(debugging) << "[DEBUG] Searching for CEMRGNET output";

    while (searchit.hasNext()) {
        QFileInfo searchfinfo(searchit.next());
        if (searchfinfo.fileName().contains(".nii", Qt::CaseSensitive)) {
            if (searchfinfo.fileName().contains("dcm-LGE", Qt::CaseSensitive))
                lgePath = searchfinfo.absoluteFilePath();

            if (searchfinfo.fileName().contains("dcm-MRA", Qt::CaseSensitive))
                mraPath = searchfinfo.absoluteFilePath();

            if (debugging && searchfinfo.fileName().contains("LA-cemrgnet.nii", Qt::CaseSensitive))
                cnnPath = searchfinfo.absoluteFilePath();
        }
    }//_while

    QDialog* inputs = new QDialog(0, 0);
    m_UIcemrgnet.setupUi(inputs);
    connect(m_UIcemrgnet.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIcemrgnet.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    QString meType_UI;
    int minStep_UI = -1, maxStep_UI = 3;
    int methodType_UI = 2, thresh_methodType_UI = 1;
    QStringList separated_thresh_list;
    std::vector<double> values_vector;

    if (dialogCode == QDialog::Accepted) {

        MITK_INFO << "[UI] User inputs being selected.";
        MITK_INFO << "[UI] Intensity projection";
        bool ok1, ok2;
        minStep_UI = m_UIcemrgnet.minStep_lineEdit->text().toInt(&ok1);
        maxStep_UI = m_UIcemrgnet.maxStep_lineEdit->text().toInt(&ok2);
        if (!ok1) minStep_UI = -1;
        if (!ok2) maxStep_UI = 3;

        methodType_UI = m_UIcemrgnet.maxProjection_radioButton->isChecked() ? 2 : 1;
        meType_UI = m_UIcemrgnet.maxProjection_radioButton->isChecked() ? "Max" : "Mean";

        MITK_INFO << ("[UI] Using: " + meType_UI + " Intensity projection.").toStdString();
        MITK_INFO << ("[UI] In/out values: (" + QString::number(minStep_UI) + ", " +
            QString::number(maxStep_UI) + ")").toStdString();
        MITK_INFO << QString::number(methodType_UI);
        MITK_INFO << "[UI] Thresholding information.";

        //bool ok3;
        thresh_methodType_UI = m_UIcemrgnet.iir_radioButton->isChecked() ? 1 : 2;
        QString thresh_list, whichThresh;

        if (m_UIcemrgnet.iir_radioButton->isChecked()) { // IIR method
            whichThresh = "IIR";
            thresh_list = m_UIcemrgnet.iir_textEdit->toPlainText();
            separated_thresh_list << "0.97" << "1.16";
        } else if (m_UIcemrgnet.meanSD_radioButton->isChecked()) { // SDev method
            whichThresh = "MEAN+SD";
            thresh_list = m_UIcemrgnet.meanSD_textEdit->toPlainText();
            separated_thresh_list << "2.3" << "3.3";
        }

        MITK_INFO << ("[UI] Threshold: " + whichThresh).toStdString();
        MITK_INFO << ("[UI] Threshold list: " + thresh_list).toStdString();
        MITK_INFO << QString::number(thresh_methodType_UI);

        thresh_list.remove(" ", Qt::CaseSensitive);
        if (!thresh_list.isEmpty()) {

            MITK_INFO << "[UI] Creating list of thresholds";
            separated_thresh_list.removeLast();
            separated_thresh_list.removeLast();
            separated_thresh_list = thresh_list.split(",", QString::SkipEmptyParts);
            int listspaces = separated_thresh_list.removeAll(" ");
            int listduplicates = separated_thresh_list.removeDuplicates();
            separated_thresh_list.sort();
            if (debugging) {
                MITK_INFO << ("[UI][DEBUG] Spaces: " + QString::number(listspaces)).toStdString();
                MITK_INFO << ("[UI][DEBUG] Duplicates: " + QString::number(listduplicates)).toStdString();
            }

        }//_if

        for (int ix = 0; ix < separated_thresh_list.size(); ix++) {
            MITK_INFO << separated_thresh_list.at(ix);
            bool vOK;
            double tryNumber = separated_thresh_list.at(ix).toDouble(&vOK);
            if (vOK) values_vector.push_back(tryNumber);
        }

        inputs->close();
        inputs->deleteLater();

    } else if (dialogCode == QDialog::Rejected) {

        MITK_INFO << "[ATTENTION] Cancelled automatic analysis.";
        QMessageBox::warning(
            NULL, "Automatic analysis cancelled",
            "'Cancel' button pressed, no calculations were made.");
        inputs->close();
        inputs->deleteLater();
        return;

    }//_if

    MITK_INFO << ("Files to be read: \n\n [LGE]: " + lgePath + "\n [MRA]: " + mraPath).toStdString();

    if (!mraPath.isEmpty()) {

        vtkSmartPointer<vtkTimerLog> timerLog = vtkSmartPointer<vtkTimerLog>::New();
        typedef itk::Image<short, 3> ImageTypeSHRT;
        typedef itk::Image<short, 3> ImageTypeCHAR;
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        MITK_INFO << "[AUTOMATIC_ANALYSIS] Setting Docker on MIRTK to OFF";
        cmd->SetUseDockerContainers(_useDockerInPlugin);

        timerLog->StartTimer();
        if (cnnPath.isEmpty()) {
            MITK_INFO << "[AUTOMATIC_ANALYSIS] Computing automatic segmentation step.";
            cnnPath = cmd->DockerCemrgNetPrediction(mraPath);
        }

        MITK_INFO << "Round pixel values from automatic segmentation.";
        CemrgCommonUtils::RoundPixelValues(cnnPath);

        if (!cnnPath.isEmpty()) {

            MITK_INFO << ("Successful prediction with file " + cnnPath).toStdString();
            // QString direct = finfo.absolutePath();
            MITK_INFO << "[AUTOMATIC_ANALYSIS][1] Adjust CNN label to MRA";
            mitk::Image::Pointer mraIMG = mitk::IOUtil::Load<mitk::Image>(mraPath.toStdString());
            mitk::Image::Pointer cnnIMG = mitk::IOUtil::Load<mitk::Image>(cnnPath.toStdString());
            double origin[3]; double spacing[3];
            mraIMG->GetGeometry()->GetOrigin().ToArray(origin);
            mraIMG->GetGeometry()->GetSpacing().ToArray(spacing);

            vtkSmartPointer<vtkImageResize> resizeFilter = vtkSmartPointer<vtkImageResize>::New();
            resizeFilter->SetResizeMethodToOutputDimensions();
            resizeFilter->SetOutputDimensions(mraIMG->GetDimension(0), mraIMG->GetDimension(1), mraIMG->GetDimension(2));
            resizeFilter->InterpolateOff();
            resizeFilter->SetInputData(cnnIMG->GetVtkImageData());
            resizeFilter->Update();

            vtkSmartPointer<vtkImageChangeInformation> changeFilter = vtkSmartPointer<vtkImageChangeInformation>::New();
            changeFilter->SetInputConnection(resizeFilter->GetOutputPort());
            changeFilter->SetOutputSpacing(spacing);
            changeFilter->SetOutputOrigin(origin);
            changeFilter->Update();

            cnnIMG->Initialize(changeFilter->GetOutput());
            cnnIMG->SetVolume(changeFilter->GetOutput()->GetScalarPointer());

            MITK_INFO << "[AUTOMATIC_ANALYSIS][2] Image registration";
            cnnPath = direct + "/LA.nii";
            QString laregPath = direct + "/LA-reg.nii";

            mitk::IOUtil::Save(cnnIMG, cnnPath.toStdString());
            cmd->ExecuteRegistration(direct, lgePath, mraPath); // rigid.dof is the default name
            cmd->ExecuteTransformation(direct, cnnPath, laregPath);

            MITK_INFO << "[AUTOMATIC_ANALYSIS][3] Clean segmentation";
            typedef itk::ImageRegionIteratorWithIndex<ImageTypeCHAR> ItType;
            typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
            typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
            using DuplicatorType = itk::ImageDuplicator<ImageTypeCHAR>;

            ImageTypeCHAR::Pointer orgSegImage = ImageTypeCHAR::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(laregPath.toStdString()), orgSegImage);

            ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
            connected1->SetInput(orgSegImage);
            connected1->Update();

            LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr1 = LabelShapeKeepNObjImgFilterType::New();
            lblShpKpNObjImgFltr1->SetInput(connected1->GetOutput());
            lblShpKpNObjImgFltr1->SetBackgroundValue(0);
            lblShpKpNObjImgFltr1->SetNumberOfObjects(1);
            lblShpKpNObjImgFltr1->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
            lblShpKpNObjImgFltr1->Update();

            DuplicatorType::Pointer duplicator = DuplicatorType::New();
            duplicator->SetInputImage(lblShpKpNObjImgFltr1->GetOutput());
            duplicator->Update();
            ItType itDUP(duplicator->GetOutput(), duplicator->GetOutput()->GetRequestedRegion());
            for (itDUP.GoToBegin(); !itDUP.IsAtEnd(); ++itDUP)
                if ((int)itDUP.Get() != 0)
                    itDUP.Set(1);
            QString segCleanPath = direct + "/prodClean.nii";
            mitk::IOUtil::Save(mitk::ImportItkImage(duplicator->GetOutput()), segCleanPath.toStdString());
            MITK_INFO << ("[...][3.1] Saved file: " + segCleanPath).toStdString();

            MITK_INFO << "[AUTOMATIC_ANALYSIS][4] Vein clipping mesh";
            QString output1 = cmd->ExecuteSurf(direct, segCleanPath, "close", 1, .5, 0, 10);
            mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(output1.toStdString());
            vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
            deci->SetInputData(shell->GetVtkPolyData());
            deci->SetTargetReduction(0.1);
            deci->PreserveTopologyOn();
            deci->Update();
            shell->SetVtkPolyData(deci->GetOutput());

            vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
            vtkSmartPointer<vtkPolyData> pd = shell->Clone()->GetVtkPolyData();
            for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
                double* point = pd->GetPoint(i);
                point[0] = -point[0];
                point[1] = -point[1];
                pd->GetPoints()->SetPoint(i, point);
            }//_for
            pointLocator->SetDataSet(pd);
            pointLocator->BuildLocator();

            MITK_INFO << "[AUTOMATIC_ANALYSIS][5] Separate veins";
            typedef itk::BinaryCrossStructuringElement<ImageTypeCHAR::PixelType, 3> CrossType;
            typedef itk::BinaryMorphologicalOpeningImageFilter<ImageTypeCHAR, ImageTypeCHAR, CrossType> MorphFilterType;
            typedef itk::RelabelComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> RelabelFilterType;

            ImageTypeCHAR::Pointer veinsSegImage = lblShpKpNObjImgFltr1->GetOutput();
            ItType itORG(orgSegImage, orgSegImage->GetRequestedRegion());
            ItType itVEN(veinsSegImage, veinsSegImage->GetRequestedRegion());
            itORG.GoToBegin();

            for (itVEN.GoToBegin(); !itVEN.IsAtEnd(); ++itVEN) {
                if ((int)itVEN.Get() != 0)
                    itVEN.Set((int)itORG.Get());
                ++itORG;
            }
            for (itVEN.GoToBegin(); !itVEN.IsAtEnd(); ++itVEN)
                if ((int)itVEN.Get() != 2)
                    itVEN.Set(0);

            CrossType binaryCross;
            binaryCross.SetRadius(2.0);
            binaryCross.CreateStructuringElement();
            MorphFilterType::Pointer morphFilter = MorphFilterType::New();
            morphFilter->SetInput(veinsSegImage);
            morphFilter->SetKernel(binaryCross);
            morphFilter->SetForegroundValue(2);
            morphFilter->SetBackgroundValue(0);
            morphFilter->UpdateLargestPossibleRegion();
            veinsSegImage = morphFilter->GetOutput();
            mitk::IOUtil::Save(mitk::ImportItkImage(veinsSegImage), (direct + "/prodVeins.nii").toStdString());

            ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();
            connected2->SetInput(veinsSegImage);
            connected2->Update();

            RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
            relabeler->SetInput(connected2->GetOutput());
            relabeler->Update();
            mitk::IOUtil::Save(mitk::ImportItkImage(relabeler->GetOutput()), (direct + "/prodSeparatedVeins.nii").toStdString());
            MITK_INFO << ("[...][5.1] Saved file: " + direct + "/prodSeparatedVeins.nii").toStdString();

            MITK_INFO << "[AUTOMATIC_ANALYSIS][6] Find vein landmark";
            veinsSegImage = relabeler->GetOutput();
            ItType itLMK(veinsSegImage, veinsSegImage->GetRequestedRegion());
            vtkSmartPointer<vtkIdList> pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
            pickedSeedIds->Initialize();
            std::vector<std::vector<double>> veinsCentre;
            const int nveins = static_cast<int>(connected2->GetObjectCount());

            MITK_INFO << ("[...][6.1] Number of veins found: " + QString::number(nveins)).toStdString();
            for (int j = 0; j < nveins; j++) {
                int ctrVeinsVoxels = 0;
                std::vector<double> veinLandmark(3, 0.0);
                for (itLMK.GoToBegin(); !itLMK.IsAtEnd(); ++itLMK) {
                    if ((int)itLMK.Get() == (j + 1)) {
                        ImageTypeCHAR::PointType point;
                        veinsSegImage->TransformIndexToPhysicalPoint(itLMK.GetIndex(), point);
                        veinLandmark[0] += point[0];
                        veinLandmark[1] += point[1];
                        veinLandmark[2] += point[2];
                        ctrVeinsVoxels++;
                    }
                }//_for
                veinLandmark[0] /= ctrVeinsVoxels;
                veinLandmark[1] /= ctrVeinsVoxels;
                veinLandmark[2] /= ctrVeinsVoxels;
                veinsCentre.push_back(veinLandmark);
            }//_nveins
            for (int j = 0; j < nveins; j++) {
                double veinLandmark[3];
                veinLandmark[0] = veinsCentre.at(j)[0];
                veinLandmark[1] = veinsCentre.at(j)[1];
                veinLandmark[2] = veinsCentre.at(j)[2];
                vtkIdType id = pointLocator->FindClosestPoint(veinLandmark);
                pickedSeedIds->InsertNextId(id);
            }//_nveins
            std::vector<int> pickedSeedLabels;
            for (int j = 0; j < nveins; j++)
                pickedSeedLabels.push_back(21);

            MITK_INFO << "[AUTOMATIC_ANALYSIS][7] Clip the veins";

            std::unique_ptr<CemrgAtriaClipper> clipper(new CemrgAtriaClipper(direct, shell));
            bool successful = clipper->ComputeCtrLines(pickedSeedLabels, pickedSeedIds, true);
            if (!successful) {
                QMessageBox::critical(NULL, "Attention", "Computation of Centrelines Failed!");
                return;
            }//_Check for failure
            MITK_INFO << "[...][7.1] ComputeCtrLines finished .";

            successful = clipper->ComputeCtrLinesClippers(pickedSeedLabels);
            if (!successful) {
                QMessageBox::critical(NULL, "Attention", "Computation of Clipper Planes Failed!");
                return;
            }//_if
            MITK_INFO << "[...][7.2] ComputeCtrLinesClippers finished .";

            clipper->ClipVeinsImage(pickedSeedLabels, mitk::ImportItkImage(duplicator->GetOutput()), false);
            MITK_INFO << "[...][7.3] ClipVeinsImage finished .";

            MITK_INFO << "[AUTOMATIC_ANALYSIS][8] Create a mesh from clipped segmentation of veins";
            QString output2 = cmd->ExecuteSurf(direct, (direct + "/PVeinsCroppedImage.nii"), "close", 1, .5, 0, 10);
            mitk::Surface::Pointer LAShell = mitk::IOUtil::Load<mitk::Surface>(output2.toStdString());

            MITK_INFO << "[AUTOMATIC_ANALYSIS][9] Clip the mitral valve";
            ImageTypeCHAR::Pointer mvImage = ImageTypeCHAR::New();
            mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(segCleanPath.toStdString()), mvImage);
            ItType itMVI1(mvImage, mvImage->GetRequestedRegion());
            itORG.GoToBegin();
            for (itMVI1.GoToBegin(); !itMVI1.IsAtEnd(); ++itMVI1) {
                if ((int)itMVI1.Get() != 0)
                    itMVI1.Set((int)itORG.Get());
                ++itORG;
            }
            for (itMVI1.GoToBegin(); !itMVI1.IsAtEnd(); ++itMVI1)
                if ((int)itMVI1.Get() != 3)
                    itMVI1.Set(0);
            typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
            ConnectedComponentImageFilterType::Pointer connected3 = ConnectedComponentImageFilterType::New();
            connected3->SetInput(mvImage);
            connected3->Update();
            typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
            LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr2 = LabelShapeKeepNObjImgFilterType::New();
            lblShpKpNObjImgFltr2->SetInput(connected3->GetOutput());
            lblShpKpNObjImgFltr2->SetBackgroundValue(0);
            lblShpKpNObjImgFltr2->SetNumberOfObjects(1);
            lblShpKpNObjImgFltr2->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
            lblShpKpNObjImgFltr2->Update();
            mvImage = lblShpKpNObjImgFltr2->GetOutput();
            mitk::IOUtil::Save(mitk::ImportItkImage(mvImage), (direct + "/prodMVI.nii").toStdString());

            // Make vtk of prodMVI
            QString mviShellPath = cmd->ExecuteSurf(direct, "prodMVI.nii", "close", 1, 0.5, 0, 10);
            // Implement code from command line tool
            mitk::Surface::Pointer ClipperSurface = mitk::IOUtil::Load<mitk::Surface>(mviShellPath.toStdString());
            vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFn = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
            implicitFn->SetInput(ClipperSurface->GetVtkPolyData());
            vtkMTimeType mtime = implicitFn->GetMTime();
            std::cout << "MTime: " << mtime << std::endl;
            vtkSmartPointer<vtkClipPolyData> mvclipper = vtkSmartPointer<vtkClipPolyData>::New();
            mvclipper->SetClipFunction(implicitFn);
            mvclipper->SetInputData(LAShell->GetVtkPolyData());
            mvclipper->InsideOutOff();
            mvclipper->Update();

            MITK_INFO << "[...][9.1] Extract and clean surface mesh.";
            vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSe