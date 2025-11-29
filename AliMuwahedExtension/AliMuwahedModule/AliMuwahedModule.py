import logging
import os

import vtk

import slicer, vtk, qt, SampleData
from slicer.ScriptedLoadableModule import *
from slicer.util import *

#%%
#
# AliMuwahed
#

class AliMuwahedModule(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "AliMuwahed"  # TODO: make this more human readable by adding spaces
        self.parent.categories = ["Examples"]  # TODO: set categories (folders where the module shows up in the module selector)
        self.parent.dependencies = []  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["Caroline Essert (University of Strasbourg)"]  # TODO: replace with "Firstname Lastname (Organization)"
        # TODO: update with short description of the module and a link to online module documentation
        self.parent.helpText = """
            This is an example of scripted loadable module bundled in an extension.
            See more information in <a href="https://github.com/organization/projectname#AliMuwahed">module documentation</a>.
            """
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = """
            This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
            and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
            """

#%%
#
# AliMuwahedWidget
#

class AliMuwahedWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)  # needed for parameter node observation
        # store logic in a member variable
        self.logic = AliMuwahedLogic()


    def setup(self):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        # initialisation that needs to be here - don't remove
        ScriptedLoadableModuleWidget.setup(self)

        # HelloWorld button
        helloWorldButton = qt.QPushButton("Hello world")
        helloWorldButton.toolTip = "Print 'Hello world' in standard window"
        self.layout.addWidget(helloWorldButton)
        helloWorldButton.connect('clicked(bool)', self.onHelloWorldButtonClicked)

        # PrintPos button
        PrintPosButton = qt.QPushButton("Print pos")
        self.layout.addWidget(PrintPosButton)
        PrintPosButton.connect('clicked(bool)', self.onPrintPosButtonButtonClicked)


    # HelloWorld button callback function
    def onHelloWorldButtonClicked(self):
        # get message to display from the logic
        message = self.logic.process()
        # display the message in a separate window (message box)
        qt.QMessageBox.information(slicer.util.mainWindow(), 'Slicer Python', message)

    # PrintPos button callback function
    def onPrintPosButtonButtonClicked(self):
        # call logic function to print position of the first fiducial (requires a fiducial to be added beforehand to work properly)
        self.logic.printPosF1()

#%%
#
# AliMuwahedLogic
#

class AliMuwahedLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)

    def process(self):
        return "Hello world!"

    def printPosF1(self, caller=None, event=None): 
        # gets node that contains all fiducials (usually called "F")
        try:
            f = getNode('F')
            # initilize a position
            pos=[0,0,0]
            # get the first fiducial's position in pos (fiducial of index=0)
            f.GetNthControlPointPosition(0,pos)
            # print the position coordinates in the python console
            print(pos)
        except slicer.util.MRMLNodeNotFoundException:
            print("Please create a fiducial first")


#%%
#
# AliMuwahedTest
#

class AliMuwahedTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """ Do whatever is needed to reset the state - typically a scene clear will be enough.
        """
        # open the MRBrainTumor1 sample dataset
        sampleDataLogic = SampleData.SampleDataLogic()
        masterVolumeNode = sampleDataLogic.downloadMRBrainTumor1()

    def runTest(self):
        """Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_AliMuwahed1()

    def test_AliMuwahed1(self):
        """ Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        # quick message box to inform that the test is starting
        self.delayDisplay("Starting the test")

        # create node to store fiducials
        markupsNodeID = slicer.modules.markups.logic().AddNewFiducialNode()
        markupsNode = getNode(markupsNodeID)
        markupsNode.SetName("F")
        # add one fiducial
        markupsNode.AddFiducial(6.4, 35.1, 0.7)

        # get the logic
        logic = AliMuwahedLogic()
        # call function printPosF1 to test it on the previously added fiducial
        logic.printPosF1()

        # quick message box to inform that the test has successfully ended
        self.delayDisplay('Test passed')
