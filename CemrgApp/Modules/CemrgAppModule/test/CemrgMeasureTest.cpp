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
    QTest::newRow("3-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}} << -1.0;
    QTest::newRow("Test 1") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}} << 1.7320508075688772935274463415059;
    QTest::newRow("Test 2") << Points{{7.0, 4.0, 3.0}, {17.0, 6.0, 2.0}} << 10.246950765959598383221038680521;
    QTest::newRow("Test 3") << Points{{-4.0, 3.0, -2.0}, {23.0, -2.0, 5.0}} << 28.337254630609507934884031143657;
}

void TestCemrgMeasure::CalcDistance() {
    QFETCH(Points, points);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->CalcDistance(points), result);
}

void TestCemrgMeasure::CalcPerimeter_data() {
    QTest::addColumn<Points>("points");
    QTest::addColumn<double>("result");

    QTest::newRow("1-point test") << Points{{0.0, 0.0, 0.0}} << -1.0;
    QTest::newRow("2-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}} << -1.0;
    QTest::newRow("3-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}} << 3.4641016151377545870548926830117;
    QTest::newRow("4-point test") << Points{{-1.0, 5.0, 7.0}, {2.0, 3.0, 4.0}, {21.0, -15.0, -20.0}, {-1.0, 5.0, 7.0}} << 80.363148824999240938914956776348;
    QTest::newRow("5-point test") << Points{{-13.0, 8.0, -7.0}, {1.0, 1.0, -1.0}, {-21.0, 15.0, 20.0}, {12.0, 51.0, 72.0}, {-13.0, 8.0, -7.0}} << 214.93578434989057588768898364186;
}

void TestCemrgMeasure::CalcPerimeter() {
    QFETCH(Points, points);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->CalcPerimeter(points), result);
}

void TestCemrgMeasure::CalcArea_data() {
    QTest::addColumn<Points>("points");
    QTest::addColumn<double>("result");

    QTest::newRow("1-point test") << Points{{0.0, 0.0, 0.0}} << -1.0;
    QTest::newRow("2-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}} << -1.0;
    QTest::newRow("3-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}} << 0.0000000160123397681644855524243;
    QTest::newRow("4-point test") << Points{{-1.0, 5.0, 7.0}, {2.0, 3.0, 4.0}, {21.0, -15.0, -20.0}, {-1.0, 5.0, 7.0}} << 11.3688170009044213770721398759633;
    QTest::newRow("5-point test") << Points{{-13.0, 8.0, -7.0}, {1.0, 1.0, -1.0}, {-21.0, 15.0, 20.0}, {12.0, 51.0, 72.0}, {-13.0, 8.0, -7.0}} << 1016.3855905904007386197918094694614;
}

void TestCemrgMeasure::CalcArea() {
    QFETCH(Points, points);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->CalcArea(points), result);
}

void TestCemrgMeasure::FindCentre_data() {
    QTest::addColumn<mitk::PointSet::Pointer>("pointSet");
    QTest::addColumn<mitk::Point3D>("result");

    // Point and Result values
    const array<tuple<Point, Point>, 8> findCentreData { {
        { {0, 0, 0}, {0, 0, 0} },
        { {1, 1, 1}, {0.5, 0.5, 0.5} },
        { {-1, 5, 7}, {0, 2, 2.66666666666666651864} },
        { {2, 3, 4}, {0.5, 2.25, 3} },
        { {21, -15, -20}, {4.6, -1.2, -1.6} },
        { {-13, 8, -7}, {1.66666666666666674068, 0.33333333333333331483, -2.5} },
        { {12, 51, 72}, {3.14285714285714279370, 7.57142857142857117481, 8.1