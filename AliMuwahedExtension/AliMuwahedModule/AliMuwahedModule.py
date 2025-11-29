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
        # Q7: Ablation metrics display section (collapsible)
        self.ablationCollapsible = slicer.qMRMLCollapsibleButton()
        self.ablationCollapsible.text = "Ablation Metrics"
        self.layout.addWidget(self.ablationCollapsible)
        self.ablationGrid = qt.QGridLayout()
        self.ablationCollapsible.setLayout(self.ablationGrid)
        self.ablationLabels = []  # Store references to ablation metric labels

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
        # Q7: Update ablation metrics after autoplace
        self.logic.updateAblationMetrics(self)
        
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
        # Slicer Markups API does not provide GetSelectedControlPoint().
        # Instead, we can use GetNumberOfControlPoints() and always show the last needle moved.
        n = caller.GetNumberOfControlPoints()
        if n < 2:
            return
        # Show distances for the last needle (most recently added/moved)
        needleIdx = (n - 1) // 2
        self.logic.computeSingleNeedleVesselDistances(self, needleIdx)

    def onAutoPlaceButtonClicked(self):
        # Q2: Automatic needle placement
        # Moves the tip of each needle to the center of mass of the tumor mesh.
        self.logic.autoPlaceNeedleTip(self)
        # Q7: Update ablation metrics after autoplace
        self.logic.updateAblationMetrics(self)

#%%
#
# AliMuwahedModuleLogic
#

class AliMuwahedModuleLogic(ScriptedLoadableModuleLogic):
    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.needleLineSources = []  # Store line sources for each needle
        self.needleModels = []       # Store model nodes for each needle
        self.ablationSpheres = []    # Store ablation sphere model nodes for each needle
        self.fiducialNode = None

    # --- Mesh/Labelmap helpers (Q7) ---
    def modelNodeToLabelmap(self, modelNode, referenceVolumeNode=None, spacing=1.0):
        """
        Converts a model node (tumor or ablation sphere) to a binary labelmap (vtkImageData).
        If referenceVolumeNode is provided, uses its geometry; otherwise, creates a new image covering the model bounds.
        Spacing is in mm.
        Returns vtkImageData (binary: 1 inside, 0 outside).
        """
        polyData = modelNode.GetPolyData()
        bounds = polyData.GetBounds()
        minX, maxX, minY, maxY, minZ, maxZ = bounds
        dimX = int((maxX - minX) / spacing) + 1
        dimY = int((maxY - minY) / spacing) + 1
        dimZ = int((maxZ - minZ) / spacing) + 1
        image = vtk.vtkImageData()
        image.SetSpacing(spacing, spacing, spacing)
        image.SetOrigin(minX, minY, minZ)
        image.SetDimensions(dimX, dimY, dimZ)
        image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
        image.GetPointData().GetScalars().Fill(0)
        pol2stenc = vtk.vtkPolyDataToImageStencil()
        pol2stenc.SetInputData(polyData)
        pol2stenc.SetOutputOrigin(image.GetOrigin())
        pol2stenc.SetOutputSpacing(image.GetSpacing())
        pol2stenc.SetOutputWholeExtent(image.GetExtent())
        pol2stenc.Update()
        imgstenc = vtk.vtkImageStencil()
        imgstenc.SetInputData(image)
        imgstenc.SetStencilData(pol2stenc.GetOutput())
        imgstenc.ReverseStencilOff()
        imgstenc.SetBackgroundValue(0)
        imgstenc.Update()
        outImage = imgstenc.GetOutput()
        outImage.GetPointData().GetScalars().SetName("Label")
        return outImage

    def countVoxelsInLabelmap(self, labelmap, labelValue=1):
        """
        Counts the number of voxels in a vtkImageData labelmap that have the given labelValue (default 1).
        Returns the count and the volume in mm^3 (count * voxel volume).
        """
        arr = labelmap.GetPointData().GetScalars()
        count = 0
        for i in range(arr.GetNumberOfTuples()):
            if arr.GetValue(i) == labelValue:
                count += 1
        spacing = labelmap.GetSpacing()
        voxelVolume = spacing[0] * spacing[1] * spacing[2]
        totalVolume = count * voxelVolume
        return count, totalVolume

    # --- Geometry/Update logic ---
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
        # Q6: Create ablation sphere at 5mm from internal tip
        ablationRadius = 10.0  # 10mm diameter (example, adjust as needed)
        direction = [pt2[j] - pt1[j] for j in range(3)]
        length = sum([d**2 for d in direction]) ** 0.5
        if length == 0:
            offset = [0,0,5]
        else:
            offset = [direction[j]/length*5.0 for j in range(3)]
        sphereCenter = [pt1[j] + offset[j] for j in range(3)]
        sphereSource = vtk.vtkSphereSource()
        sphereSource.SetCenter(*sphereCenter)
        sphereSource.SetRadius(ablationRadius/2.0)
        sphereSource.SetThetaResolution(32)
        sphereSource.SetPhiResolution(32)
        sphereSource.Update()
        sphereModelNode = slicer.modules.models.logic().AddModel(sphereSource.GetOutput())
        sphereModelNode.SetName(f"Ablation-{len(self.ablationSpheres)+1}")
        sphereModelNode.GetDisplayNode().SetColor(1,0,0)  # Red
        sphereModelNode.GetDisplayNode().SetOpacity(0.3)
        self.ablationSpheres.append(sphereModelNode)
        fiducialNode.AddObserver(vtk.vtkCommand.ModifiedEvent, lambda caller, event: self.onFiducialMoved(caller, event, widget))
        fiducialNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, lambda caller, event: self.onFiducialMoved(caller, event, widget))

    def updateAllNeedlesFromFiducials(self, caller, event):
        # Q1/Q4/Q6: Update needle geometry and ablation sphere interactively
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
            direction = [pt2[j] - pt1[j] for j in range(3)]
            length = sum([d**2 for d in direction]) ** 0.5
            if length == 0:
                offset = [0,0,5]
            else:
                offset = [direction[j]/length*5.0 for j in range(3)]
            sphereCenter = [pt1[j] + offset[j] for j in range(3)]
            ablationRadius = 10.0
            sphereSource = vtk.vtkSphereSource()
            sphereSource.SetCenter(*sphereCenter)
            sphereSource.SetRadius(ablationRadius/2.0)
            sphereSource.SetThetaResolution(32)
            sphereSource.SetPhiResolution(32)
            sphereSource.Update()
            self.ablationSpheres[i].SetAndObservePolyData(sphereSource.GetOutput())

    def onFiducialMoved(self, caller, event, widget):
        # Q4: Automatic update of the distance
        self.updateAllNeedlesFromFiducials(caller, event)
        self.computeNeedleVesselDistances(widget)
        # Q7: Update ablation metrics live when fiducial moves
        self.updateAblationMetrics(widget)

    # --- Tumor/Needle logic ---
    def autoPlaceNeedleTip(self, widget):
        # Q2: Automatic needle placement
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

    # --- Risk/Distance logic ---
    def computeNeedleVesselDistances(self, widget):
        # Q3/Q5: Computation of distances and automatic coloring according to risk
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
            if minDistForNeedle < minRisk:
                minRisk = minDistForNeedle
                minRiskNeedleIdx = needleIdx
        for i, modelNode in enumerate(self.needleModels):
            if i == minRiskNeedleIdx:
                modelNode.GetDisplayNode().SetColor(1,0,0)
            else:
                modelNode.GetDisplayNode().SetColor(1,1,0)
        for label in getattr(widget, 'distanceLabels', []):
            widget.distanceGrid.removeWidget(label)
            label.deleteLater()
        widget.distanceLabels = []
        header = qt.QLabel("Needle/Vessel")
        widget.distanceGrid.addWidget(header, 0, 0)
        widget.distanceLabels.append(header)
        for v, vessel_name in enumerate(vessel_names):
            vesselLabel = qt.QLabel(vessel_name)
            widget.distanceGrid.addWidget(vesselLabel, 0, v+1)
            widget.distanceLabels.append(vesselLabel)
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
        for label in getattr(widget, 'distanceLabels', []):
            widget.distanceGrid.removeWidget(label)
            label.deleteLater()
        widget.distanceLabels = []
        header = qt.QLabel(f"Needle {needleIdx+1}")
        widget.distanceGrid.addWidget(header, 0, 0)
        widget.distanceLabels.append(header)
        for v, vessel_name in enumerate(vessel_names):
            vesselLabel = qt.QLabel(vessel_name)
            widget.distanceGrid.addWidget(vesselLabel, 0, v+1)
            widget.distanceLabels.append(vesselLabel)
        for v, value in enumerate(needle_results):
            distLabel = qt.QLabel(value)
            widget.distanceGrid.addWidget(distLabel, 1, v+1)
            widget.distanceLabels.append(distLabel)

    # --- Q7: Ablation metrics ---
    def computeAblationTumorOverlap(self, tumorNode, ablationSphereNodes, spacing=1.0):
        """
        Computes overlap between tumor and ablation spheres.
        Returns: ablated tumor voxel count, ablated tumor volume, total tumor voxel count, total tumor volume, efficiency (fraction).
        Steps:
        1. Convert tumor mesh to binary labelmap (voxels inside tumor = 1).
        2. Convert each ablation sphere mesh to binary labelmap.
        3. Union all ablation sphere labelmaps (any voxel covered by any sphere = 1).
        4. For each voxel, if inside both tumor and ablation, count as ablated.
        5. Count all tumor voxels (labelmap = 1).
        6. Compute volumes by multiplying voxel count by voxel volume (spacing³).
        7. Efficiency = ablated tumor volume / total tumor volume.
        """
        # Step 1: Convert tumor mesh to labelmap
        tumorLabelmap = self.modelNodeToLabelmap(tumorNode, spacing=spacing)
        # Step 2 & 3: Convert ablation spheres to labelmaps and union them
        ablationLabelmap = None
        for sphereNode in ablationSphereNodes:
            sphereLabelmap = self.modelNodeToLabelmap(sphereNode, spacing=spacing)
            if ablationLabelmap is None:
                ablationLabelmap = vtk.vtkImageData()
                ablationLabelmap.DeepCopy(sphereLabelmap)
            else:
                # Union: set voxels to 1 if either is 1
                arrA = ablationLabelmap.GetPointData().GetScalars()
                arrB = sphereLabelmap.GetPointData().GetScalars()
                for i in range(arrA.GetNumberOfTuples()):
                    if arrB.GetValue(i) == 1:
                        arrA.SetValue(i, 1)
        # Step 4: Intersection - tumor voxels inside ablation
        arrTumor = tumorLabelmap.GetPointData().GetScalars()
        arrAblation = ablationLabelmap.GetPointData().GetScalars()
        ablatedCount = 0
        for i in range(arrTumor.GetNumberOfTuples()):
            if arrTumor.GetValue(i) == 1 and arrAblation.GetValue(i) == 1:
                ablatedCount += 1
        # Step 5: Count all tumor voxels
        totalCount, totalVolume = self.countVoxelsInLabelmap(tumorLabelmap, labelValue=1)
        # Step 6: Compute ablated tumor volume
        spacing = tumorLabelmap.GetSpacing()
        voxelVolume = spacing[0] * spacing[1] * spacing[2]
        ablatedVolume = ablatedCount * voxelVolume
        # Step 7: Compute efficiency
        efficiency = ablatedVolume / totalVolume if totalVolume > 0 else 0
        return ablatedCount, ablatedVolume, totalCount, totalVolume, efficiency

    def updateAblationMetrics(self, widget):
        """
        Updates the ablation metrics UI section with current ablation/tumor overlap and efficiency.
        Steps:
        1. Get tumor and ablation sphere nodes.
        2. Compute overlap and metrics using computeAblationTumorOverlap.
        3. Clear previous UI labels.
        4. Display metrics: ablated voxels, ablated volume, total tumor voxels, total tumor volume, efficiency.
        """
        try:
            tumorNode = getNode('livertumor04')
        except slicer.util.MRMLNodeNotFoundException:
            return
        ablationSphereNodes = self.ablationSpheres
        if not ablationSphereNodes:
            return
        # Step 2: Compute metrics
        ablatedCount, ablatedVolume, totalCount, totalVolume, efficiency = self.computeAblationTumorOverlap(tumorNode, ablationSphereNodes, spacing=1.0)
        # Step 3: Clear previous labels
        for label in getattr(widget, 'ablationLabels', []):
            widget.ablationGrid.removeWidget(label)
            label.deleteLater()
        widget.ablationLabels = []
        # Step 4: Add metrics to UI
        metrics = [
            ("Ablated Tumor Voxels", ablatedCount),
            ("Ablated Tumor Volume (mm³)", f"{ablatedVolume:.2f}"),
            ("Total Tumor Voxels", totalCount),
            ("Total Tumor Volume (mm³)", f"{totalVolume:.2f}"),
            ("Ablation Efficiency (%)", f"{efficiency*100:.2f}")
        ]
        for i, (labelText, value) in enumerate(metrics):
            label = qt.QLabel(f"{labelText}: {value}")
            widget.ablationGrid.addWidget(label, i, 0)
            widget.ablationLabels.append(label)
        # Optionally, print to console for debug
        print(f"Ablation metrics: Ablated {ablatedCount} voxels, {ablatedVolume:.2f} mm³; Total {totalCount} voxels, {totalVolume:.2f} mm³; Efficiency {efficiency*100:.2f}%")
