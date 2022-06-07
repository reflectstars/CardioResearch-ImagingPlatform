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

#include "CemrgStrainsTest.hpp"

static bool FuzzyCompare(const double& lhs, const double& rhs) {
    constexpr double almostZero = 1e-6;
    if (abs(lhs) < almostZero && abs(rhs) < almostZero)
        return true;
    
    return qFuzzyCompare(lhs, rhs);
}

mitk::DataNode::Pointer TestCemrgStrains::ReferenceAHA(const array<int, 3>& segRatios, bool pacingSite) {
    mitk::PointSet::Pointer pointSet = mitk::IOUtil::Load<mitk::PointSet>((QFINDTESTDATA(CemrgTestData::strainPath) + "/PointSet.mps").toStdString());
    mitk::DataNode::Pointer lmNode = mitk::DataNode::New();
    lmNode->SetData(pointSet);
    cemrgStrains->ReferenceAHA(lmNode, (int*)segRatios.data(), pacingSite);
    return lmNode;
}

void TestCemrgStrains::initTestCase() {

}

void TestCemrgStrains::cleanupTestCase() {

}

void TestCemrgStrains::CalculateGlobalSqzPlot_data() {
    QTest::addColumn<int>("meshNo");
    QTest::addColumn<double>("result");

    const array<double, CemrgTestData::strainDataSize> globalSqzPlotData {
        0,
        -0.08891526089125015297
    };

    for (size_t i = 0; i < globalSqzPlotData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << (int)i << globalSqzPlotData[i];
}

void TestCemrgStrains::CalculateGlobalSqzPlot() {
    QFETCH(int, meshNo);
    QFETCH(double, result);

    QCOMPARE(cemrgStrains->CalculateGlobalSqzPlot(meshNo), result);
}

void TestCemrgStrains::CalculateSqzPlot_data() {
    QTest::addColumn<int>("meshNo");
    QTest::addColumn<vector<double>>("result");

    const array<vector<double>, CemrgTestData::strainDataSize> sqzPlotData { {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {-0.24161484861895057841, -0.24203718408602067913, -0.24189488048215221361, -0.24269039612066256595, -0.23629833084106455221, -0.23554401643925146348, -0.23956905910100148582, -0.24574671778317599968, -0.24545802466742852599, -0.24045512800893234506, -0.23522178212071434555, -0.23486487750507906158, -0.24153056562196148493, -0.23671392457799511622, -0.24116868663725185562, -0.23596466822878597869}
    } };

    // Preparation for tests
    ReferenceAHA();

    for (size_t i = 0; i < sqzPlotData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << (int)i << sqzPlotData[i];
}

void TestCemrgStrains::CalculateSqzPlot() {
    QFETCH(int, meshNo);
    QFETCH(vector<double>, result);

    QVERIFY(equal(begin(result), end(result), begin(cemrgStrains->CalculateSqzPlot(meshNo)), FuzzyCompare));
}

void TestCemrgStrains::CalculateStrainsPlot_data() {
    QTest::addColumn<int>("meshNo");
    QTest::addColumn<mitk::DataNode::Pointer>("lmNode");
    QTest::addColumn<int>("flag");
    QTest::addColumn<vector<double>>("result");

    const array<vector<double>, CemrgTestData::strainDataSize> strainsPlotData { {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {-0.45874886139987286482, -0.38998985665445839999, -0.39095134301528899901, -0.45061430609189162544, -0.48968384803220571522, -0.49392625117662009027, -0.45850249287033251200, -0.35280369788968640732, -0.35498114762612859030, -0.44960627315573070684, -0.49382093507974794688, -0.49665076176258193819, -0.31543513441107573492, -0.43