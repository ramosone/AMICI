#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupSteadystate)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupSteadystate, testSimulation) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/nosensi/");
    simulateAndWriteToFile(model, "/model_steadystate/nosensi/");
    delete model;
}

TEST(groupSteadystate, testSensitivityForward) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/sensiforward/");
    simulateAndWriteToFile(model, "/model_steadystate/sensiforward/");
    delete model;
}

TEST(groupSteadystate, testSensitivityForwardPlist) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/sensiforwardplist/");
    simulateAndWriteToFile(model, "/model_steadystate/sensiforwardplist/");
    delete model;
}


TEST(groupSteadystate, testSensitivityForwardErrorInt) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/sensiforwarderrorint/");
    simulateAndWriteToFile(model, "/model_steadystate/sensiforwarderrorint/");
    delete model;
}

TEST(groupSteadystate, testSensitivityForwardErrorNewt) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/sensiforwarderrornewt/");
    simulateAndWriteToFile(model, "/model_steadystate/sensiforwarderrornewt/");
    delete model;
}


TEST(groupSteadystate, testSensitivityForwardDense) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/sensiforwarddense/");
    delete model;
}

TEST(groupSteadystate, testSensitivityForwardSPBCG) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_steadystate/nosensiSPBCG/");
    delete model;
}


