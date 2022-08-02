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
          