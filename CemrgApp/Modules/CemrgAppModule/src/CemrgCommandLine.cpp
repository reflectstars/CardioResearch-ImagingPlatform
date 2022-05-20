
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
 * Commandline Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>

// Qt
#include <QFileDialog>
#include <QFileInfo>
#include <QFile>
#include <QCoreApplication>
#include <QDebug>
#include <QDir>
#include <QMessageBox>

// C++ Standard
#include <thread>
#include <chrono>
#include <sys/stat.h>
#include "CemrgCommandLine.h"

CemrgCommandLine::CemrgCommandLine() {

    _useDockerContainers = true;
    _debugvar = false;
    _dockerimage = "biomedia/mirtk:v1.1.0";

    //Setup panel
    panel = new QTextEdit(0,0);
    QPalette palette = panel->palette();
    palette.setColor(QPalette::Base, Qt::black);
    palette.setColor(QPalette::Text, Qt::red);
    panel->setPalette(palette);
    panel->setReadOnly(true);

    //Setup dialog
    layout = new QVBoxLayout();
    dial = new QDialog(0,0);
    dial->setFixedSize(640, 480);
    dial->setLayout(layout);
    dial->layout()->addWidget(panel);
    dial->show();

    //Setup the process
    process = std::unique_ptr<QProcess>(new QProcess(this));
    process->setProcessChannelMode(QProcess::MergedChannels);
    connect(process.get(), SIGNAL(readyReadStandardOutput()), this, SLOT(UpdateStdText()));
    connect(process.get(), SIGNAL(readyReadStandardError()), this, SLOT(UpdateErrText()));
    connect(process.get(), SIGNAL(finished(int)), this, SLOT(FinishedAlert()));
}

CemrgCommandLine::~CemrgCommandLine() {

    process->close();
    dial->deleteLater();
    panel->deleteLater();
    layout->deleteLater();
}

QDialog* CemrgCommandLine::GetDialog() {

    return dial;
}

/***************************************************************************
 ****************** Execute Plugin Specific Functions **********************
 ***************************************************************************/

QString CemrgCommandLine::ExecuteSurf(QString dir, QString segPath, QString morphOperation, int iter, float th, int blur, int smth) {

    MITK_INFO << "[ATTENTION] SURFACE CREATION: Close -> Surface -> Smooth";

    QString closeOutputPath, surfOutputPath;
    QString outAbsolutePath = "ERROR_IN_PROCESSING";
    closeOutputPath = ExecuteMorphologicalOperation(morphOperation, dir, segPath, "segmentation.s.nii", iter);

    mitk::ProgressBar::GetInstance()->Progress();
    if (QString::compare(closeOutputPath, "ERROR_IN_PROCESSING")!=0) {

        surfOutputPath = ExecuteExtractSurface(dir, closeOutputPath, "segmentation.vtk", th, blur);
        mitk::ProgressBar::GetInstance()->Progress();

        if (QString::compare(surfOutputPath, "ERROR_IN_PROCESSING")!=0) {

            outAbsolutePath = ExecuteSmoothSurface(dir, surfOutputPath, surfOutputPath, smth);
            remove((dir + "/segmentation.s.nii").toStdString().c_str());
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            mitk::ProgressBar::GetInstance()->Progress();
        }//_if
    } else {
        mitk::ProgressBar::GetInstance()->Progress(2);
    }//_if

    return outAbsolutePath;
}

QString CemrgCommandLine::ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName) {

    MITK_INFO << "[ATTENTION] Attempting MeshTools3D libraries.";

    QString segmentationDirectory = dir + "/";
    QString outputDirectory = segmentationDirectory + "CGALMeshDir";
    QString outAbsolutePath = outputDirectory + "/" + outputName + ".vtk"; // many outputs are created with meshtools3d. .vtk is the one used in CemrgApp

    MITK_INFO << "Using static MeshTools3D libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/M3DLib";
    QString executableName = executablePath + "/meshtools3d";
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << "-f" << paramsFullPath;
        arguments << "-seg_dir" << segmentationDirectory;;
        arguments << "-seg_name" << segmentationName;
        arguments << "-out_dir" << outputDirectory;
        arguments << "-out_name" << outputName;

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MeshTools3D libraries not found");
        MITK_WARN << "MeshTools3D libraries not found. Please make sure the M3DLib folder is inside the directory:\n\t" + mitk::IOUtil::GetProgramPath();
    }//_if

    //Setup EnVariable - in windows TBB_NUM_THREADS should be set in the system environment variables
#ifndef _WIN32
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("TBB_NUM_THREADS","12");
    process->setProcessEnvironment(env);
#endif

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful) {
        MITK_WARN << "MeshTools3D did not produce a good outcome.";
        return "ERROR_IN_PROCESSING";
    } else {
        return outAbsolutePath;
    }
}

void CemrgCommandLine::ExecuteTracking(QString dir, QString imgTimes, QString param, QString output) {

    MITK_INFO << "[ATTENTION] Attempting Registration.";

    QString commandName = "register";
    QString imgTimesFilePath, outAbsolutePath;
    QString prodPath = dir + "/";

    imgTimesFilePath = imgTimes.contains(dir, Qt::CaseSensitive) ? imgTimes : prodPath + imgTimes;
    outAbsolutePath = output.contains(dir, Qt::CaseSensitive) ? output : prodPath + output;

    if (!outAbsolutePath.contains(".dof", Qt::CaseSensitive)) outAbsolutePath += ".dof";

    MITK_INFO << ("[...] IMAGE FILES PATH: " + imgTimesFilePath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << "-images" << imgTimesFilePath;
        if (!param.isEmpty()) arguments << "-parin" << param;
        arguments << "-dofout" << outAbsolutePath;
        arguments << "-threads" << "12";
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                        mitk::IOUtil::GetProgramPath();
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {

    //Time Smoothness
    int fctTime = 10;
    noFrames *= smooth;
    if (smooth == 2) {
        fctTime = 5;
    } else if (smooth == 5) {
        fctTime = 2;
    }

    QString output = dir + "/transformed-";
    QString thisOutput;
    for (int i=0; i<noFrames; i++) {
        thisOutput = output + QString::number(i)+".vtk";
        ExecuteTransformationOnPoints(dir, inputMesh, thisOutput, dofin, iniTime);
        iniTime += fctTime;
        mitk::ProgressBar::GetInstance()->Progress();
    }
}

void CemrgCommandLine::ExecuteRegistration(QString dir, QString fixed, QString moving, QString transformFileName, QString modelname) {

    MITK_INFO << "[ATTENTION] Attempting Registration.";

    //lge : fixed   ||   mra : moving
    QString commandName = "register";
    QString fixedfullpath, movingfullpath, outAbsolutePath;
    QString prodPath = dir + "/";

    fixedfullpath = fixed.contains(dir, Qt::CaseSensitive) ? fixed : prodPath + fixed;
    movingfullpath = moving.contains(dir, Qt::CaseSensitive) ? moving : prodPath + moving;
    outAbsolutePath = transformFileName.contains(dir, Qt::CaseSensitive) ? transformFileName : prodPath + transformFileName;

    if (!fixedfullpath.contains(".nii", Qt::CaseSensitive)) fixedfullpath += ".nii";
    if (!movingfullpath.contains(".nii", Qt::CaseSensitive)) movingfullpath += ".nii";
    if (!outAbsolutePath.contains(".dof", Qt::CaseSensitive)) movingfullpath += ".dof";

    MITK_INFO << ("[...] MOVING (source): " + movingfullpath).toStdString();
    MITK_INFO << ("[...] FIXED (target): " + fixedfullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << movingfullpath;
        arguments << fixedfullpath;
        arguments << "-dofout" << outAbsolutePath;
        arguments << "-model" << modelname;
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+mitk::IOUtil::GetProgramPath();
    }//_if

    MITK_INFO << ("Performing a " + modelname + " registration").toStdString();
    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteTransformation(QString dir, QString imgname, QString regname, QString transformFileFullPath) {

    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString commandName = "transform-image";
    QString dofpath, imgNamefullpath, outAbsolutePath;
    QString prodPath = dir + "/";

    imgNamefullpath = imgname.contains(dir, Qt::CaseSensitive) ? imgname : prodPath + imgname;
    outAbsolutePath = regname.contains(dir, Qt::CaseSensitive) ? regname : prodPath + regname;
    dofpath = transformFileFullPath.contains(dir, Qt::CaseSensitive) ? transformFileFullPath : prodPath + transformFileFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << imgNamefullpath; //input
        arguments << outAbsolutePath; //output
        arguments << "-dofin" << dofpath;
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+mitk::IOUtil::GetProgramPath();
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString transformFileName, bool transformThePoints) {

    MITK_INFO << "[ATTENTION] Attempting INIT-DOF.";

    QString executablePath, executableName, commandName, sourceMeshPath, targetMeshPath, outAbsolutePath;
    QString prodPath = dir + "/";
    QStringList arguments;

    commandName = "init-dof"; //simple translation
    sourceMeshPath = sourceMeshP.contains(dir, Qt::CaseSensitive) ? sourceMeshP : prodPath + sourceMeshP;
    targetMeshPath = targetMeshP.contains(dir, Qt::CaseSensitive) ? targetMeshP : prodPath + targetMeshP;
    outAbsolutePath = transformFileName.contains(dir, Qt::CaseSensitive) ? transformFileName : prodPath + transformFileName;

    MITK_INFO << "Using static MIRTK libraries.";
    executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);

    if (apathd.exists()) {
