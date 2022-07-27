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
    bool test = std::adjacent_find(indexNodes.begin(), indexNodes.