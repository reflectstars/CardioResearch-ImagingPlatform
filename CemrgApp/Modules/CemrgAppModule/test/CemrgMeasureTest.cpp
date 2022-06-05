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
 * CEMRGAPPMODULE TESTS
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#include "CemrgMeasureTest.hpp"

typedef CemrgMeasure::Point Point;
typedef CemrgMeasure::Points Points;

void TestCemrgMeasure::initTestCase() {
    // Create surface data
    for (size_t i = 0; i < surfaceData.size(); i++) {
        surfaceData[i].first = QFINDTESTDATA(CemrgTestData::surfacePaths[i]);
        surfaceData[i].second = mitk::IOUtil::Load<mitk::Surface>(surfaceData[i].first.toStdString());
    }
}

void TestCemrgMeasure::cleanupTestCase() {

}

void TestCemrgMeasure::CalcDistance_data() {
    QTest::addColumn<Points>("points");
    QTest::addColumn<double>("result");

    QTest::newRow("1-point test") << Points{{0.0, 0.0, 0.0}} << -1.0;
    QTest::newRow