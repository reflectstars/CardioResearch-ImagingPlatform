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
#include <QDebug>
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
    parser.setDescription(
        "This template command-line app for testing different functionalities.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    parser.addArgument(
        "reference", "i", mitkCommandLineParser::InputFile,
        "Input object (reference)", "Reference for registration. Accepts any image format known to MITK.",
        us::Any(), false);
    parser.addArgument(
        "otherimage", "j", mitkCommandLineParser::InputFile,
        "Input object (reference)", "Image to calculate transformation from. Accepts any image format known to MITK.",
        us::Any(), false);
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["reference"].Empty() ||
        parsedArgs["otherimage"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto input1 = us::any_cast<std::string>(parsedArgs["reference"]);
    auto input2 = us::any_cast<std::string>(parsedArgs["otherimage"]);

    // Default values for optional arguments
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose"))
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);

    try {
        // Code the functionality of the cmd app here.
        if (verbose) 