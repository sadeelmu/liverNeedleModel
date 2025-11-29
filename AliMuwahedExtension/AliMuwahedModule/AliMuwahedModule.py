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

        # Observe selection changes on fiducials (use only ModifiedEvent)
        try:
            fiducialNode = getNode('F')
            fiducialNode.AddObserver(vtk.vtkCommand.ModifiedEvent, self.onFiducialSelected)
        except slicer.util.MRMLNodeNotFoundException:
            pass

    # Create Needles button callback function
    def onCreateNeedlesButtonClicked(self):
        # Q1: Add new needle (pair of fiducials and cylinder)
        self.logic.createNeedles(self)

    # Automatic Needle Placement button
    def onAutoPlaceButtonClicked(self):
        # Q2: Place needle tips at tumor center
        self.logic.autoPlaceNeedleTip(self)
        
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
    def createNeedles(self, widget):
        # Get or create fiducial node
        try:
            fiducialNode = getNode('F')
        except slicer.util.MRMLNodeNotFoundException:
            fiducialNode = slicer.modules.markups.logic().AddNewFiducialNode()
            fiducialNode = getNode(fiducialNode)
            fiducialNode.SetName('F')
        self.fiducialNode = fiducialNode
        fiducialNode.GetDisplayNode().SetGlyphScale(3.0)
        n = fiducialNode.GetNumberOfControlPoints()
        pt1 = [0, 0, 0]
        pt2 = [0, 0, 150]
        fiducialNode.AddControlPoint(vtk.vtkVector3d(*pt1))
        fiducialNode.SetNthControlPointLabel(n, f"F-{n+1}")
        fiducialNode.AddControlPoint(vtk.vtkVector3d(*pt2))
        fiducialNode.SetNthControlPointLabel(n+1, f"F-{n+2}")
        # VTK pipeline: LineSource -> TubeFilter -> TriangleFilter
        lineSource = vtk.vtkLineSource()
        lineSource.SetPoint1(pt1)
        lineSource.SetPoint2(pt2)
        tubeFilter = vtk.vtkTubeFilter()
        tubeFilter.SetInputConnection(lineSource.GetOutputPort())
        tubeFilter.SetRadius(1.0)
        tubeFilter.SetNumberOfSides(20)
        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(tubeFilter.GetOutputPort())
        triangleFilter.Update()
        modelNode = slicer.modules.models.logic().AddModel(triangleFilter.GetOutput())
        modelNode.SetName(f"Needle-{len(self.needleLineSources)+1}")
        modelNode.GetDisplayNode().SetColor(1,1,0)
        modelNode.GetDisplayNode().SetEdgeVisibility(True)
        modelNode.GetDisplayNode().SetRepresentation(1)
        self.needleLineSources.append(lineSource)
        self.needleModels.append(modelNode)
        # Pass widget to observer using lambda for live update
        fiducialNode.AddObserver(vtk.vtkCommand.ModifiedEvent, lambda caller, event: self.onFiducialMoved(caller, event, widget))
        # Add real-time observer for point movement
        fiducialNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, lambda caller, event: self.onFiducialMoved(caller, event, widget))

    def onFiducialMoved(self, caller, event, widget):
        # Update all needle geometries interactively when any fiducial point is moved
        self.updateAllNeedlesFromFiducials(caller, event)
        # Automatically update distances in the UI
        self.computeNeedleVesselDistances(widget)

    def updateAllNeedlesFromFiducials(self, caller, event):
        if not self.fiducialNode:
            return
        for i, lineSource in enumerate(self.needleLineSources):
            idx1, idx2 = i*2, i*2+1
            if self.fiducialNode.GetNumberOfControlPoints() <= idx2:
                continue
            pt1 = [0,0,0]
            pt2 = [0,0,0]
            self.fiducialNode.GetNthControlPointPosition(idx1, pt1)
            self.fiducialNode.GetNthControlPointPosition(idx2, pt2)
            lineSource.SetPoint1(pt1)
            lineSource.SetPoint2(pt2)
            tubeFilter = vtk.vtkTubeFilter()
            tubeFilter.SetInputConnection(lineSource.GetOutputPort())
            tubeFilter.SetRadius(1.0)
            tubeFilter.SetNumberOfSides(20)
            triangleFilter = vtk.vtkTriangleFilter()
            triangleFilter.SetInputConnection(tubeFilter.GetOutputPort())
            triangleFilter.Update()
            self.needleModels[i].SetAndObservePolyData(triangleFilter.GetOutput())

    # Q2: Automatic needle tip placement at tumor center of mass
    def autoPlaceNeedleTip(self, widget):
        try:
            tumorNode = getNode('livertumor04')
            polyData = tumorNode.GetPolyData()
            centerOfMassFilter = vtk.vtkCenterOfMass()
            centerOfMassFilter.SetInputData(polyData)
            centerOfMassFilter.Update()
            com = centerOfMassFilter.GetCenter()
        except:
            return
        if self.fiducialNode:
            for i in range(0, self.fiducialNode.GetNumberOfControlPoints(), 2):
                if self.fiducialNode.GetNumberOfControlPoints() > i+1:
                    pt2 = [0,0,0]
                    self.fiducialNode.GetNthControlPointPosition(i+1, pt2)
                    old_pt1 = [0,0,0]
                    self.fiducialNode.GetNthControlPointPosition(i, old_pt1)
                    direction = [pt2[j] - old_pt1[j] for j in range(3)]
                    self.fiducialNode.SetNthControlPointPosition(i, com)
                    new_pt2 = [com[j] + direction[j] for j in range(3)]
                    self.fiducialNode.SetNthControlPointPosition(i+1, new_pt2)
        self.updateAllNeedlesFromFiducials(self.fiducialNode, None)
        self.computeNeedleVesselDistances(widget)

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
        Also update needle colors: closest needle to any vessel is colored red, others yellow.
        """
        vessel_names = [
            'portalvein', 'venoussystem', 'artery'
        ]
        results = []
        minRisk = float('inf')
        minRiskNeedleIdx = None
        for needleIdx, modelNode in enumerate(self.needleModels):
            needlePolyData = modelNode.GetPolyData()
            needle_results = []
            minDistForNeedle = float('inf')
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
                if minDist < minDistForNeedle:
                    minDistForNeedle = minDist
            results.append(needle_results)
            # Track needle with smallest risk (closest to any vessel)
            if minDistForNeedle < minRisk:
                minRisk = minDistForNeedle
                minRiskNeedleIdx = needleIdx
        # Update needle colors
        for i, modelNode in enumerate(self.needleModels):
            if i == minRiskNeedleIdx:
                modelNode.GetDisplayNode().SetColor(1,0,0)  # Red
            else:
                modelNode.GetDisplayNode().SetColor(1,1,0)  # Yellow
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
