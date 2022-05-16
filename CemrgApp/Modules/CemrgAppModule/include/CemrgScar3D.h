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
 * Scar Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgScar3D_h
#define CemrgScar3D_h

// Qmitk
#include <mitkImage.h>
#include <mitkSurface.h>
#include <mitkPointSet.h>
#include <vtkFloatArray.h>
#include <MitkCemrgAppModuleExports.h>
#include <QString>

class MITKCEMRGAPPMODULE_EXPORT CemrgScar3D {

public:

    CemrgScar3D();
    mitk::Surface::Pointer Sc