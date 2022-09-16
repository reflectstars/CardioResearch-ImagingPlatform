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

// Qmitk
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include "WallThicknessCalculationsClipperView.h"
#include "WallThicknessCalculationsView.h"

// VTK
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>

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
#include <QDesktopWidget>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>

QString WallThicknessCalculationsClipperView::fileName;
QString WallThicknessCalculationsClipperView::directory;
const std::string WallThicknessCalculationsClipperView::VIEW_ID = "org.mitk.views.wathcaclipperview";

WallThicknessCalculationsClipperView::WallThicknessCalculationsClipperView() {
    this->inputs = new QDialog(0, 0);
}

void WallThicknessCalculationsClipperView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(CtrLines()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(CtrPlanes()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(ClipperImage()));
    connect(m_Controls.slider, SIGNAL(valueChanged(int)), this, SLOT(CtrPlanesPlacer()));
    connect(m_Controls.spinBox, SIGNAL(valueChanged(double)), this, SLOT(CtrPlanesPlacer()));
    connect(m_Controls.comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(CtrLinesSelector(int)));

    //Create GUI widgets
    inputs = new QDialog(0, 0);
    m_Labels.setupUi(inputs);
    connect(m_Labels.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_Labels.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

    //Setup renderer
    surfActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0, 0, 0);
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    std::string shortcuts = "R: reset centrelines\nSpace: add seed point\nDelete: remove seed point";
    txtActor->SetInput(shortcuts.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    m_Controls.widget_1->SetRenderWindow(renderWindow);
    m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

    //Setup keyboard interactor
    callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    callBack->SetCallback(KeyCallBackFunc);
    callBack->SetClientData(this);
    interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    interactor->GetInteractorStyle()->KeyPressActivationOff();
    interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);

    //Initialisation
    iniPreSurf();
    if (surface.IsNotNull()) {
        pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
        pickedSeedIds->Initialize();
        pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
        pickedLineSeeds->Initialize();
        pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
        pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
        pickedCutterSeeds->Initialize();
        pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
        clipper = std::unique_ptr<CemrgAtriaClipper>(new CemrgAtriaClipper(directory, surface));
        Visualiser();
    }
    m_Controls.slider->setEnabled(false);
    m_Controls.spinBox->setEnabled(false);
    m_Controls.comboBox->setEnabled(false);
}

void WallThicknessCalculationsClipperView::SetFocus() {
    m_Controls.button_1->setFocus();
}

void WallThicknessCalculationsClipperView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

WallThicknessCalculationsClipperView::~WallThicknessCalculationsClipperView() {
    inputs->deleteLater();
}

void WallThicknessCalculationsClipperView::SetDirectoryFile(const QString directory, const QString fileName) {
    WallThicknessCalculationsClipperView::fileName = fileName;
    WallThicknessCalculationsClipperView::directory = directory;
}

void WallThicknessCalculationsClipperView::iniPreSurf() {
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to clip!");
        this->GetSite()->GetPage()->ResetPerspective();
        return;
    }

    //Find the selected node
    QString path;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            //Check seg node name
            if (segNode->GetName().compare(fileName.left(fileName.lastIndexOf(QChar('.'))).toStdString()) != 0) {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                this->GetSite()->GetPage()->ResetPerspective();
                return;
            }//_if

            //Ask for user input to set the parameters
            QDialog* inputs = new QDialog(0, 0);
            m_UIMeshing.setupUi(inputs);
            connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
            connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
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
                    QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
                if (!ok1) th = 0.5;
                if (!ok2) bl = 0.8;
                if (!ok3) smth = 3;
                if (!ok4) ds = 0.5;
                //_if

                //Mesh creation
                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(2);
                auto filter = mitk::ManualSegmentationToSurfaceFilter::New();
                filter->SetInput(image);
                filter->SetThreshold(th);
                filter->SetUseGaussianImageSmooth(true);
                filter->SetSmooth(true);
                filter->SetMedianFilter3D(true);
                filter->InterpolationOn();
                filter->SetGaussianStandardDeviation(bl);
                filter->SetMedianKernelSize(smth, smth, smth);
                filter->SetDecimate(mitk::ImageToSurfaceFilter::QuadricDecimation);
                filter->SetTargetReduction(ds);
                filter->UpdateLargestPossibleRegion();
                mitk::ProgressBar::GetInstance()->Progress();
                mitk::Surface::Pointer shell = filter->GetOutput();
                vtkSmartPointer<vtkPolyData> pd = shell->GetVtkPolyData();
                pd->SetVerts(nullptr);
                pd->SetLines(nullptr);
                for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
                    double* point = pd->GetPoint(i);
                    point[0] = -point[0];
                    point[1] = -point[1];
                    pd->GetPoints()->SetPoint(i, point);
                }//_for
                vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
                connectivityFilter->SetInputData(pd);
                connectivityFilter->ColorRegionsOff();
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
                normals->AutoOrientNormalsOn();
                normals->FlipNormalsOff();
                normals->SetInputConnection(connectivityFilter->GetOutputPort());
                normals->Update();
                shell->SetVtkPolyData(normals->GetOutput());
                surfac