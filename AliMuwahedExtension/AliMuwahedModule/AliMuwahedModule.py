import logging
import os

import vtk

import slicer, vtk, qt, SampleData
from slicer.ScriptedLoadableModule import *
from slicer.util import *

#%%
#
# AliMuwahedModule
#

class AliMuwahedModule(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "AliMuwahedModule"  # TODO: make this more human readable by adding spaces
        self.parent.categories = ["Examples"]  # TODO: set categories (folders where the module shows up in the module selector)
        self.parent.dependencies = []  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["Caroline Essert (University of Strasbourg)"]  # TODO: replace with "Firstname Lastname (Organization)"
        # TODO: update with short description of the module and a link to online module documentation
        self.parent.helpText = """
            This is an example of scripted loadable module bundled in an extension.
            See more information in <a href="https://github.com/organization/projectname#AliMuwahedModule">module documentation</a>.
            """
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = """
            This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
            and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
            """

#%%
#
# AliMuwahedModuleWidget
#

class AliMuwahedModuleWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
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
        self.logic = AliMuwahedModuleLogic()

    def setup(self):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        # initialisation that needs to be here - don't remove
        ScriptedLoadableModuleWidget.setup(self)

        # PrintPos button
        PrintPosButton = qt.QPushButton("Print pos")
        self.layout.addWidget(PrintPosButton)
        PrintPosButton.connect('clicked(bool)', self.onPrintPosButtonButtonClicked)

        # Create Needles button: Adds a new needle (pair of fiducials and cylinder) each time
        createNeedlesButton = qt.QPushButton("Create Needles")
        createNeedlesButton.toolTip = "Add a new pair of fiducials and cylinder (needle) each time you press."
        self.layout.addWidget(createNeedlesButton)
        createNeedlesButton.connect('clicked(bool)', self.onCreateNeedlesButtonClicked)

        # Q2: Automatic Needle Placement button
        autoPlaceButton = qt.QPushButton("Auto Place Needle Tip")
        autoPlaceButton.toolTip = "Automatically place needle tip (F-1, F-3, ...) at tumor center of mass."
        self.layout.addWidget(autoPlaceButton)
        autoPlaceButton.connect('clicked(bool)', self.onAutoPlaceButtonClicked)

        # Q3: Compute Distances button
        computeDistancesButton = qt.QPushButton("Compute Needle-Vessel Distances")
        computeDistancesButton.toolTip = "Compute and display minimum distance from each needle to each vessel."
        self.layout.addWidget(computeDistancesButton)
        computeDistancesButton.connect('clicked(bool)', self.onComputeDistancesButtonClicked)

        # Q3: Distance display section (collapsible)
        self.distanceCollapsible = slicer.qMRMLCollapsibleButton()
        self.distanceCollapsible.text = "Distance"
        self.layout.addWidget(self.distanceCollapsible)
        self.distanceGrid = qt.QGridLayout()
        self.distanceCollapsible.setLayout(self.distanceGrid)
        self.distanceLabels = []  # Store references to labels for updating

        # Observe selection changes on fiducials
        try:
            fiducialNode = getNode('F')
            fiducialNode.AddObserver(vtk.vtkCommand.ControlPointModifiedEvent, self.onFiducialSelected)
        except slicer.util.MRMLNodeNotFoundException:
            pass

    # Create Needles button callback function
    def onCreateNeedlesButtonClicked(self):
        # Q1: Add new needle (pair of fiducials and cylinder)
        self.logic.createNeedles()

    # Automatic Needle Placement button
    def onAutoPlaceButtonClicked(self):
        # Q2: Place needle tips at tumor center
        self.logic.autoPlaceNeedleTip()
        
    # PrintPos button callback function
    def onPrintPosButtonButtonClicked(self):
        # call logic function to print position of the first fiducial (requires a fiducial to be added beforehand to work properly)
        self.logic.printPosF1()

    def onComputeDistancesButtonClicked(self):
        # Q3: Compute and display needle-vessel distances in the interface
        self.logic.computeNeedleVesselDistances(self)

    def onFiducialSelected(self, caller, event):
        # Find selected fiducial index
        selectedIndex = caller.GetSelectionNode().GetActiveMarkupsFiducialID()
        if selectedIndex is None:
            return
        # Find which needle this fiducial belongs to
        idx = caller.GetSelectedControlPoint()  # This should give the index
        if idx is None:
            return
        needleIdx = idx // 2  # Each needle has 2 fiducials
        self.logic.computeSingleNeedleVesselDistances(self, needleIdx)

#%%
#
# AliMuwahedModuleLogic
#

class AliMuwahedModuleLogic(ScriptedLoadableModuleLogic):
    """Implements all computation for the module."""
    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.needleLineSources = []  # Store line sources for each needle
        self.needleModels = []       # Store model nodes for each needle
        self.fiducialNode = None

    # Q1: Needle creation and interaction
    def createNeedles(self):
        """
        Add a new needle (cylinder) between a new pair of fiducial points each time the button is pressed.
        Does not remove previous needles or fiducials.
        """
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

        # Determine index for new fiducials
        n = fiducialNode.GetNumberOfControlPoints()
        pt1 = [0, 0, 0]
        pt2 = [0, 0, 150]  # 15cm along Z
        fiducialNode.AddFiducial(*pt1)
        fiducialNode.SetNthFiducialLabel(n, f"F-{n+1}")
        fiducialNode.AddFiducial(*pt2)
        fiducialNode.SetNthFiducialLabel(n+1, f"F-{n+2}")

        # VTK pipeline: LineSource -> TubeFilter -> TriangleFilter
        lineSource = vtk.vtkLineSource()
        lineSource.SetPoint1(pt1)
        lineSource.SetPoint2(pt2)
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
        needleIndex = len(self.needleLineSources) + 1
        modelNode = slicer.modules.models.logic().AddModel(triangleFilter.GetOutput())
        modelNode.SetName(f"Needle-{needleIndex}")
        modelNode.GetDisplayNode().SetColor(1,1,0)  # Yellow
        modelNode.GetDisplayNode().SetEdgeVisibility(True)
        modelNode.GetDisplayNode().SetRepresentation(1)  # Wireframe

        self.needleLineSources.append(lineSource)
        self.needleModels.append(modelNode)

        # Observe fiducial movement to update all needles and distances automatically
        fiducialNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.onFiducialMoved)
        fiducialNode.AddObserver(vtk.vtkCommand.PointModifiedEvent, self.onFiducialMoved)

    def onFiducialMoved(self, caller, event):
        # Automatically update distances when any fiducial is moved
        # Only update if a 'free' point (F-2, F-4, ...) is moved
        numPoints = caller.GetNumberOfControlPoints()
        # Check which points are moved (for simplicity, update for any move)
        # You can optimize to only update for even indices if needed
        self.computeNeedleVesselDistances(self.widget)

    def updateAllNeedlesFromFiducials(self, caller, event):
        """
        Update all needle geometries interactively when any fiducial point is moved.
        """
        if not self.fiducialNode or not self.needleLineSources:
            return
        numNeedles = len(self.needleLineSources)
        for i in range(numNeedles):
            idx1 = i*2
            idx2 = idx1+1
            if self.fiducialNode.GetNumberOfControlPoints() <= idx2:
                continue
            pt1 = [0,0,0]
            pt2 = [0,0,0]
            self.fiducialNode.GetNthControlPointPosition(idx1, pt1)
            self.fiducialNode.GetNthControlPointPosition(idx2, pt2)
            lineSource = self.needleLineSources[i]
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

            self.needleModels[i].SetAndObservePolyData(triangleFilter.GetOutput())

    # Q2: Automatic needle tip placement at tumor center of mass
    def autoPlaceNeedleTip(self):
        """
        Automatically place the tip of each needle (F-1, F-3, ...) at the center of mass of the tumor mesh.
        """
        tumorNode = None
        try:
            tumorNode = getNode('livertumor04')
        except slicer.util.MRMLNodeNotFoundException:
            slicer.util.errorDisplay("Tumor model 'livertumor04' not found in scene.")
            return
        polyData = tumorNode.GetPolyData()
        if not polyData:
            slicer.util.errorDisplay("Tumor model has no polydata.")
            return
        centerOfMassFilter = vtk.vtkCenterOfMass()
        centerOfMassFilter.SetInputData(polyData)
        centerOfMassFilter.SetUseScalarsAsWeights(False)
        centerOfMassFilter.Update()
        com = centerOfMassFilter.GetCenter()
        # Move F-1, F-3, ... to center of mass
        if self.fiducialNode:
            for i in range(0, self.fiducialNode.GetNumberOfControlPoints(), 2):
                self.fiducialNode.SetNthControlPointPosition(i, com)

    def printPosF1(self, caller=None, event=None): 
        try:
            f = getNode('F')
            pos=[0,0,0]
            f.GetNthControlPointPosition(0,pos)
            print(pos)
        except slicer.util.MRMLNodeNotFoundException:
            print("Please create a fiducial first")

    def computeNeedleVesselDistances(self, widget):
        """
        Compute and display minimum distance from each needle to each vessel mesh in the module interface.
        Display results in a grid layout under a collapsible section.
        """
        vessel_names = [
            'portalvein', 'venoussystem', 'artery'
        ]
        results = []
        for needleIdx, modelNode in enumerate(self.needleModels):
            needlePolyData = modelNode.GetPolyData()
            needle_results = []
            for vessel_name in vessel_names:
                try:
                    vesselNode = getNode(vessel_name)
                except slicer.util.MRMLNodeNotFoundException:
                    needle_results.append(f"{vessel_name}: not found")
                    continue
                vesselPolyData = vesselNode.GetPolyData()
                distanceFilter = vtk.vtkDistancePolyDataFilter()
                distanceFilter.SetInputData(0, needlePolyData)
                distanceFilter.SetInputData(1, vesselPolyData)
                distanceFilter.Update()
                distances = distanceFilter.GetOutput().GetPointData().GetArray("Distance")
                if not distances:
                    needle_results.append(f"{vessel_name}: error")
                    continue
                minDist = min([distances.GetValue(i) for i in range(distances.GetNumberOfTuples())])
                needle_results.append(f"{minDist/10:.2f} cm")
            results.append(needle_results)
        # Clear previous labels
        for label in getattr(widget, 'distanceLabels', []):
            widget.distanceGrid.removeWidget(label)
            label.deleteLater()
        widget.distanceLabels = []
        # Add header
        header = qt.QLabel("Needle/Vessel")
        widget.distanceGrid.addWidget(header, 0, 0)
        widget.distanceLabels.append(header)
        for v, vessel_name in enumerate(vessel_names):
            vesselLabel = qt.QLabel(vessel_name)
            widget.distanceGrid.addWidget(vesselLabel, 0, v+1)
            widget.distanceLabels.append(vesselLabel)
        # Add data rows
        for n, needle_results in enumerate(results):
            needleLabel = qt.QLabel(f"Needle {n+1}")
            widget.distanceGrid.addWidget(needleLabel, n+1, 0)
            widget.distanceLabels.append(needleLabel)
            for v, value in enumerate(needle_results):
                distLabel = qt.QLabel(value)
                widget.distanceGrid.addWidget(distLabel, n+1, v+1)
                widget.distanceLabels.append(distLabel)

    def computeSingleNeedleVesselDistances(self, widget, needleIdx):
        """
        Compute and display minimum distance from the selected needle to each vessel mesh in the module interface.
        """
        vessel_names = [
            'portalvein', 'venoussystem', 'artery'
        ]
        if needleIdx >= len(self.needleModels):
            return
        modelNode = self.needleModels[needleIdx]
        needlePolyData = modelNode.GetPolyData()
        needle_results = []
        for vessel_name in vessel_names:
            try:
                vesselNode = getNode(vessel_name)
            except slicer.util.MRMLNodeNotFoundException:
                needle_results.append(f"{vessel_name}: not found")
                continue
            vesselPolyData = vesselNode.GetPolyData()
            distanceFilter = vtk.vtkDistancePolyDataFilter()
            distanceFilter.SetInputData(0, needlePolyData)
            distanceFilter.SetInputData(1, vesselPolyData)
            distanceFilter.Update()
            distances = distanceFilter.GetOutput().GetPointData().GetArray("Distance")
            if not distances:
                needle_results.append(f"{vessel_name}: error")
                continue
            minDist = min([distances.GetValue(i) for i in range(distances.GetNumberOfTuples())])
            needle_results.append(f"{minDist/10:.2f} cm")
        # Clear previous labels
        for label in getattr(widget, 'distanceLabels', []):
            widget.distanceGrid.removeWidget(label)
            label.deleteLater()
        widget.distanceLabels = []
        # Add header
        header = qt.QLabel(f"Needle {needleIdx+1}")
        widget.distanceGrid.addWidget(header, 0, 0)
        widget.distanceLabels.append(header)
        for v, vessel_name in enumerate(vessel_names):
            vesselLabel = qt.QLabel(vessel_name)
            widget.distanceGrid.addWidget(vesselLabel, 0, v+1)
            widget.distanceLabels.append(vesselLabel)
        # Add data row
        for v, value in enumerate(needle_results):
            distLabel = qt.QLabel(value)
            widget.distanceGrid.addWidget(distLabel, 1, v+1)
            widget.distanceLabels.append(distLabel)
