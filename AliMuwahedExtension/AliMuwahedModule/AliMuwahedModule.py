import logging
import os

import vtk

import slicer, vtk, qt, SampleData
"""
AliMuwahedModule: Slicer module for interactive thermal ablation needle planning.
Each major section is annotated with comments explaining which project question/task it implements.
"""
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
        # Q0: Module and UI setup
        # This section sets up the Slicer module widget and all UI buttons.
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
        # Q0: UI setup
        # Adds buttons for each project task and sets up the collapsible distance display grid.
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
        # Q1: Needle creation and interaction
        # Adds a new pair of fiducials and a cylinder (needle) each time the button is pressed.
        self.logic.createNeedles(self)

    # Automatic Needle Placement button
    def onAutoPlaceButtonClicked(self):
        # Q2: Automatic needle placement
        # Moves the tip of each needle to the center of mass of the tumor mesh.
        self.logic.autoPlaceNeedleTip(self)
        
    # PrintPos button callback function
    def onPrintPosButtonButtonClicked(self):
        # call logic function to print position of the first fiducial (requires a fiducial to be added beforehand to work properly)
        self.logic.printPosF1()

    def onComputeDistancesButtonClicked(self):
        # Q3: Computation of distances and display of result
        # Computes and displays minimum distance from each needle to each vessel mesh.
        self.logic.computeNeedleVesselDistances(self)

    def onFiducialSelected(self, caller, event):
        # Q3/Q4: Interactive single-needle risk display
        # When a fiducial is selected, show distances for the corresponding needle only.
        idx = caller.GetSelectedControlPoint()
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
        # Ablation sphere data (Q6)
        # Radius and offset are in the same units as the scene (mm)
        self.ablationRadius = 10.0  # default ablation radius in mm (adjustable)
        self.ablationOffsetFromTip = 5.0  # center placed 5mm from the tip
        self.ablationSphereSources = []  # vtkSphereSource for each needle
        self.ablationModels = []  # corresponding model nodes

    # Q1: Needle creation and interaction
    def createNeedles(self, widget):
        # Q1: Needle creation and interaction
        # Creates or gets the fiducial node, adds a pair of control points (needle endpoints), and builds the VTK pipeline for the cylinder.
        # Observers are added for both ModifiedEvent and PointModifiedEvent to enable live updates.
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
        # Create corresponding ablation sphere (Q6)
        sphereSource = vtk.vtkSphereSource()
        # compute initial center: 5mm from the internal tip (pt1) towards pt2
        direction = [pt2[j] - pt1[j] for j in range(3)]
        length = sum([direction[j]*direction[j] for j in range(3)]) ** 0.5
        if length > 0.0:
            center = [pt1[j] + (direction[j]/length)*self.ablationOffsetFromTip for j in range(3)]
        else:
            center = list(pt1)
        sphereSource.SetCenter(center)
        sphereSource.SetRadius(self.ablationRadius)
        sphereSource.SetThetaResolution(24)
        sphereSource.SetPhiResolution(24)
        sphereSource.Update()
        sphereModel = slicer.modules.models.logic().AddModel(sphereSource.GetOutput())
        sphereModel.SetName(f"Ablation-{len(self.ablationSphereSources)+1}")
        # semi-transparent reddish color for ablation
        sphereModel.GetDisplayNode().SetColor(1.0, 0.4, 0.4)
        sphereModel.GetDisplayNode().SetOpacity(0.35)
        self.ablationSphereSources.append(sphereSource)
        self.ablationModels.append(sphereModel)
        # Observers for live update
        fiducialNode.AddObserver(vtk.vtkCommand.ModifiedEvent, lambda caller, event: self.onFiducialMoved(caller, event, widget))
        fiducialNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, lambda caller, event: self.onFiducialMoved(caller, event, widget))

    def onFiducialMoved(self, caller, event, widget):
        # Q4: Automatic update of the distance
        # Called whenever a fiducial is moved (including during dragging). Updates needle geometry and risk display live.
        self.updateAllNeedlesFromFiducials(caller, event)
        self.computeNeedleVesselDistances(widget)

    def updateAllNeedlesFromFiducials(self, caller, event):
        # Q1/Q4: Update needle geometry interactively
        # Updates all needle geometries when any fiducial is moved.
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
            # Update corresponding ablation sphere position (Q6)
            try:
                sphereSource = self.ablationSphereSources[i]
                # compute center 5mm from internal tip (pt1) towards pt2
                direction = [pt2[j] - pt1[j] for j in range(3)]
                length = sum([direction[j]*direction[j] for j in range(3)]) ** 0.5
                if length > 0.0:
                    center = [pt1[j] + (direction[j]/length)*self.ablationOffsetFromTip for j in range(3)]
                else:
                    center = list(pt1)
                sphereSource.SetCenter(center)
                sphereSource.SetRadius(self.ablationRadius)
                sphereSource.Update()
                self.ablationModels[i].SetAndObservePolyData(sphereSource.GetOutput())
            except Exception:
                # no sphere for this needle yet
                pass

    # Q2: Automatic needle tip placement at tumor center of mass
    def autoPlaceNeedleTip(self, widget):
        # Q2: Automatic needle placement
        # Moves the tip of each needle (F-1, F-3, ...) to the center of mass of the tumor mesh, keeping direction and length.
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
        # Q3/Q5: Computation of distances and automatic coloring according to risk
        # Computes minimum distance from each needle to each vessel mesh, displays results, and colors the highest-risk needle red.
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
        # Q3/Q4: Single needle risk display
        # Computes and displays minimum distance from the selected needle to each vessel mesh in the UI.
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
