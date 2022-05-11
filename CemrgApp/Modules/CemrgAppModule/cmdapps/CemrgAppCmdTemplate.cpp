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
CEMRG CMD APP TEMPLATE
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkCommandLineParser.h>
#include <mitkIOUtil.h>

// Qt
#include <QString>
#include <QFileInfo>
#include <QProcess>

// C++ Standard
#include <algorithm>
#include <string>

// CemrgApp
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Tests");
    parser.setTitle("Template Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("This template command-line app for testing di