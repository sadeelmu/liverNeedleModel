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

#AliMuwahedModule: Metadata for the module.
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

#AliMuwahedWidget: Sets up the UI.
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
        self.layout = self.layout()  # <-- This line ensures self.layout is set

        # HelloWorld button
        helloWorldButton = qt.QPushButton("Hello world")
        helloWorldButton.toolTip = "Print 'Hello world' in standard window"
        self.layout.addWidget(helloWorldButton)
        helloWorldButton.connect('clicked(bool)', self.onHelloWorldButtonClicked)

        # PrintPos button
        PrintPosButton = qt.QPushButton("Print pos")
        self.layout.addWidget(PrintPosButton)
        PrintPosButton.connect('clicked(bool)', self.onPrintPosButtonButtonClicked)

        # Create Needles button: Automatically creates fiducials and a cylinder (needle) between them
        createNeedlesButton = qt.QPushButton("Create Needles")
        createNeedlesButton.toolTip = "Automatically create fiducials and cylinders (needles) between control points."
        self.layout.addWidget(createNeedlesButton)
        createNeedlesButton.connect('clicked(bool)', self.onCreateNeedlesButtonClicked)

    def onCreateNeedlesButtonClicked(self):
        self.logic.createNeedles()


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

#AliMuwahedLogic: Implements the actual computation.
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
        self.needleLineSources = []  # Store line sources for each needle
        self.needleModels = []       # Store model nodes for each needle
        self.fiducialNode = None
    def createNeedles(self, numNeedles=1):
        """
        Create needles (cylinders) between pairs of automatically generated fiducial points.
        - Removes previous needles and fiducials.
        - Creates two fiducials (F-1, F-2) for one needle.
        - Builds a cylinder using vtkLineSource, vtkTubeFilter, and vtkTriangleFilter.
        - Adds the cylinder to the scene in wireframe mode.
        - Increases glyph size for fiducials for better visualization.
        - Observes fiducial movement to update the needle geometry interactively.
        """
        # Remove previous needles and fiducials if any
        for modelNode in self.needleModels:
            slicer.mrmlScene.RemoveNode(modelNode)
        self.needleLineSources = []
        self.needleModels = []

        # Create or get fiducial node
        fiducialNode = None
        try:
            fiducialNode = getNode('F')
        except slicer.util.MRMLNodeNotFoundException:
            fiducialNode = slicer.modules.markups.logic().AddNewFiducialNode()
            fiducialNode = getNode(fiducialNode)
            fiducialNode.SetName('F')
        self.fiducialNode = fiducialNode

        # Increase glyph size for visibility
        fiducialNode.GetDisplayNode().SetGlyphScale(3.0)

        # Create two fiducials for one needle (can be extended for more needles)
        points = [
            [0, 0, 0],           # F-1
            [0, 0, 150]          # F-2, 15cm along Z
        ]
        fiducialNode.RemoveAllControlPoints()
        for i, pt in enumerate(points):
            fiducialNode.AddFiducial(*pt)
            fiducialNode.SetNthFiducialLabel(i, f"F-{i+1}")

        # VTK pipeline: LineSource -> TubeFilter -> TriangleFilter
        lineSource = vtk.vtkLineSource()
        lineSource.SetPoint1(points[0])
        lineSource.SetPoint2(points[1])
        lineSource.Update()

        tubeFilter = vtk.vtkTubeFilter()
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(1.0)  # 1mm radius
        tubeFilter.SetNumberOfSides(20)
        tubeFilter.Update()

        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(tubeFilter.GetOutputPort())
        triangleFilter.Update()

        # Add cylinder (needle) to scene
        modelNode = slicer.modules.models.logic().AddModel(triangleFilter.GetOutput())
        modelNode.SetName("Needle-1")
        modelNode.GetDisplayNode().SetColor(1,1,0)  # Yellow
        modelNode.GetDisplayNode().SetEdgeVisibility(True)
        modelNode.GetDisplayNode().SetRepresentation(1)  # Wireframe

        self.needleLineSources.append(lineSource)
        self.needleModels.append(modelNode)

        # Observe fiducial movement to update needle geometry
        fiducialNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.updateNeedleFromFiducials)

    def updateNeedleFromFiducials(self, caller, event):
        """
        Update the needle geometry interactively when fiducial points are moved.
        - Reads the positions of F-1 and F-2.
        - Updates the vtkLineSource, vtkTubeFilter, and vtkTriangleFilter.
        - Updates the cylinder (needle) model in the scene.
        """
        if not self.fiducialNode or not self.needleLineSources:
            return
        if self.fiducialNode.GetNumberOfControlPoints() < 2:
            return
        pt1 = [0,0,0]
        pt2 = [0,0,0]
        self.fiducialNode.GetNthControlPointPosition(0, pt1)
        self.fiducialNode.GetNthControlPointPosition(1, pt2)
        lineSource = self.needleLineSources[0]
        lineSource.SetPoint1(pt1)
        lineSource.SetPoint2(pt2)
        lineSource.Update()

        tubeFilter = vtk.vtkTubeFilter()
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(1.0)
        tubeFilter.SetNumberOfSides(20)
        tubeFilter.Update()

        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(tubeFilter.GetOutputPort())
        triangleFilter.Update()

        self.needleModels[0].SetAndObservePolyData(triangleFilter.GetOutput())

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

#AliMuwahedTest: Basic test setup, adds a fiducial and tests the print position function.
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
