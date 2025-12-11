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
        self.parent.contributors = ["Ali Kinan (University of Strasbourg), Sadel Muwahed (University of Strasbourg)"]  # TODO: replace with "Firstname Lastname (Organization)"
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
        # UI setup
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

        # Add Compute Ablation Metrics button
        computeAblationButton = qt.QPushButton("Compute Ablation Metrics")
        computeAblationButton.toolTip = "Compute ablated volume and efficiency."
        self.layout.addWidget(computeAblationButton)
        computeAblationButton.connect('clicked(bool)', self.onComputeAblationButtonClicked)
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
        # Slicer Markups API does not provide GetSelectedControlPoint().
        # Instead, we can use GetNumberOfControlPoints() and always show the last needle moved.
        n = caller.GetNumberOfControlPoints()
        if n < 2:
            return
        # Show distances for the last needle (most recently added/moved)
        needleIdx = (n - 1) // 2
        self.logic.computeSingleNeedleVesselDistances(self, needleIdx)

    def onComputeAblationButtonClicked(self):
        # Q7: Enable live ablation metrics updates
        self.logic.ablationMetricsEnabled = True
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
        self.ablationMetricsEnabled = False  # Track if ablation metrics should update live

    # --- Geometry/Update logic (Q1, Q2, Q4, Q6) ---
    def createNeedles(self, widget):
        # Q1: Needle creation and interaction
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
        ablationRadius = 10.0
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
        # Q4: Automatic update of the distance and ablation metrics
        self.updateAllNeedlesFromFiducials(caller, event)
        self.computeNeedleVesselDistances(widget)
        if self.ablationMetricsEnabled:
            self.updateAblationMetrics(widget)

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
        if self.ablationMetricsEnabled:
            self.updateAblationMetrics(widget)

    def printPosF1(self, caller=None, event=None): 
        try:
            f = getNode('F')
            pos=[0,0,0]
            f.GetNthControlPointPosition(0,pos)
            print(pos)
        except slicer.util.MRMLNodeNotFoundException:
            print("Please create a fiducial first")

    # --- Risk/Distance logic (Q3, Q5) ---
    def computeNeedleVesselDistances(self, widget):
        vessel_names = ['portalvein', 'venoussystem', 'artery']
        results = []
        
        # Compute distance to tumor for each needle to find the closest one
        minDistToTumor = float('inf')
        closestToTumorNeedleIdx = None
        try:
            tumorNode = getNode('livertumor04')
            tumorPolyData = tumorNode.GetPolyData()
            for needleIdx, modelNode in enumerate(self.needleModels):
                needlePolyData = modelNode.GetPolyData()
                distanceFilter = vtk.vtkDistancePolyDataFilter()
                distanceFilter.SetInputData(0, needlePolyData)
                distanceFilter.SetInputData(1, tumorPolyData)
                distanceFilter.Update()
                distances = distanceFilter.GetOutput().GetPointData().GetArray("Distance")
                if distances:
                    minDist = min([distances.GetValue(i) for i in range(distances.GetNumberOfTuples())])
                    if minDist < minDistToTumor:
                        minDistToTumor = minDist
                        closestToTumorNeedleIdx = needleIdx
        except slicer.util.MRMLNodeNotFoundException:
            pass
        
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
                needle_results.append(f"{minDist:.2f} mm")
            results.append(needle_results)
        
        # Color the needle closest to tumor in red, others in yellow
        for i, modelNode in enumerate(self.needleModels):
            if i == closestToTumorNeedleIdx:
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
        vessel_names = ['portalvein', 'venoussystem', 'artery']
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
            needle_results.append(f"{minDist:.2f} mm")
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

    # --- Q7: Ablation metrics (Iterative Painting) ---

    def modelNodeToLabelmap(self, inputObject, name, referenceVolumeNode=None, spacing=1.0, labelValue=1):
        """
        Robust conversion to labelmap. Handles both ModelNodes and raw PolyData.
        Ensures PADDING is added so meshes aren't clipped.
        """
        if hasattr(inputObject, "GetPolyData"):
            polyData = inputObject.GetPolyData()
        else:
            polyData = inputObject

        if polyData is None or polyData.GetNumberOfPoints() < 1:
            print(f"Error: {name} has no polyData")
            return None
            
        labelmapNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', name)
        image = vtk.vtkImageData()

        if referenceVolumeNode:
            # Case A: Align with Reference (Painting Sphere into Tumor Grid)
            refImage = referenceVolumeNode.GetImageData()
            image.DeepCopy(refImage)
            # Copy Orientation
            mat = vtk.vtkMatrix4x4()
            referenceVolumeNode.GetIJKToRASDirectionMatrix(mat)
            labelmapNode.SetIJKToRASDirectionMatrix(mat)
            labelmapNode.SetOrigin(referenceVolumeNode.GetOrigin())
            labelmapNode.SetSpacing(referenceVolumeNode.GetSpacing())
            # Initialize with 0 (Empty Canvas) to paint on
            image.GetPointData().GetScalars().Fill(0) 
        else:
            # Case B: Create New Grid (The Master Tumor Map)
            bounds = polyData.GetBounds()
            # 5mm Padding for bounding box
            padding = 5.0 
            origin = [bounds[0] - padding, bounds[2] - padding, bounds[4] - padding]
            
            dimX = int((bounds[1] - bounds[0] + (2 * padding)) / spacing) + 1
            dimY = int((bounds[3] - bounds[2] + (2 * padding)) / spacing) + 1
            dimZ = int((bounds[5] - bounds[4] + (2 * padding)) / spacing) + 1
            
            image.SetDimensions(dimX, dimY, dimZ)
            image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
            image.SetSpacing(spacing, spacing, spacing)
            image.SetOrigin(origin)
            # Initialize with 0
            image.GetPointData().GetScalars().Fill(0)

        # Create Stencil
        pol2stenc = vtk.vtkPolyDataToImageStencil()
        pol2stenc.SetInputData(polyData)
        pol2stenc.SetOutputOrigin(image.GetOrigin())
        pol2stenc.SetOutputSpacing(image.GetSpacing())
        pol2stenc.SetOutputWholeExtent(image.GetExtent())
        pol2stenc.Update()

        # Paint 1s into the image where the stencil is
        # We fill a temp image with 1s, then cut the outside to 0
        finalImage = vtk.vtkImageData()
        finalImage.DeepCopy(image)
        finalImage.GetPointData().GetScalars().Fill(labelValue) # Fill all 1
        
        imgStencil = vtk.vtkImageStencil()
        imgStencil.SetInputData(finalImage) # All 1s
        imgStencil.SetStencilData(pol2stenc.GetOutput())
        imgStencil.ReverseStencilOff()
        imgStencil.SetBackgroundValue(0) # Cut everything outside to 0
        imgStencil.Update()

        labelmapNode.SetAndObserveImageData(imgStencil.GetOutput())
        return labelmapNode

    def computeAblationTumorOverlap(self, tumorNode, ablationSphereNodes, spacing=1.0):
        # Q7: Compute volumes using Iterative Union (Robust method)
        if not tumorNode or not ablationSphereNodes:
            return 0, 0, 0, 0, 0

        # 1. Convert Tumor to Labelmap (Master Reference)
        tumorLabelmapNode = self.modelNodeToLabelmap(tumorNode, "TempTumorMap", referenceVolumeNode=None, spacing=spacing)
        if not tumorLabelmapNode:
            return 0, 0, 0, 0, 0
            
        tumorImageData = tumorLabelmapNode.GetImageData()
        
        # 2. Create "Master Ablation Image" (Accumulator)
        masterAblationImage = vtk.vtkImageData()
        masterAblationImage.DeepCopy(tumorImageData)
        masterAblationImage.GetPointData().GetScalars().Fill(0) # Start Empty
        
        # 3. Iteratively Union each sphere
        # This handles the "Union" requirement robustly
        mathFilter = vtk.vtkImageMathematics()
        mathFilter.SetOperationToMax() # 0 vs 1 -> 1

        spheresProcessed = 0
        for i, sphereNode in enumerate(ablationSphereNodes):
            # Create a labelmap for this sphere, using Tumor Grid as reference
            # This ensures they are in the exact same voxel grid
            singleSphereLabelmap = self.modelNodeToLabelmap(
                sphereNode, 
                f"TempSphereMap_{i}", 
                referenceVolumeNode=tumorLabelmapNode, 
                labelValue=1
            )
            
            if singleSphereLabelmap:
                spheresProcessed += 1
                # Union: Master = Max(Master, Sphere)
                mathFilter.SetInput1Data(masterAblationImage)
                mathFilter.SetInput2Data(singleSphereLabelmap.GetImageData())
                mathFilter.Update()
                
                masterAblationImage.DeepCopy(mathFilter.GetOutput())
                slicer.mrmlScene.RemoveNode(singleSphereLabelmap)
        
        if spheresProcessed == 0:
            slicer.mrmlScene.RemoveNode(tumorLabelmapNode)
            return 0, 0, 0, 0, 0

        # 4. Count Voxels
        tumorScalars = tumorImageData.GetPointData().GetScalars()
        ablationScalars = masterAblationImage.GetPointData().GetScalars()
        numTuples = tumorScalars.GetNumberOfTuples()

        tumorVoxelCount = 0
        ablationVoxelCount = 0
        intersectionVoxelCount = 0

        for i in range(numTuples):
            tVal = tumorScalars.GetValue(i)
            aVal = ablationScalars.GetValue(i)

            if tVal > 0:
                tumorVoxelCount += 1
            if aVal > 0:
                ablationVoxelCount += 1
            if tVal > 0 and aVal > 0:
                intersectionVoxelCount += 1

        # 5. Calculate Volumes
        sp = tumorLabelmapNode.GetSpacing()
        oneVoxelVol = sp[0] * sp[1] * sp[2]

        totalTumorVol = tumorVoxelCount * oneVoxelVol
        totalAblationVol = ablationVoxelCount * oneVoxelVol
        ablatedTumorVol = intersectionVoxelCount * oneVoxelVol

        # 6. Efficiencies
        coverage = (ablatedTumorVol / totalTumorVol) if totalTumorVol > 0 else 0
        precision = (ablatedTumorVol / totalAblationVol) if totalAblationVol > 0 else 0

        # Cleanup
        slicer.mrmlScene.RemoveNode(tumorLabelmapNode)

        return ablatedTumorVol, totalTumorVol, totalAblationVol, coverage, precision

    def updateAblationMetrics(self, widget):
        try:
            tumorNode = getNode('livertumor04')
        except slicer.util.MRMLNodeNotFoundException:
            return
            
        ablationSphereNodes = self.ablationSpheres
        if not ablationSphereNodes:
            return

        # Compute metrics
        ablatedTumorVol, totalTumorVol, totalAblationVol, coverage, precision = \
            self.computeAblationTumorOverlap(tumorNode, ablationSphereNodes, spacing=1.0)

        # Update UI
        for label in getattr(widget, 'ablationLabels', []):
            widget.ablationGrid.removeWidget(label)
            label.deleteLater()
        widget.ablationLabels = []

        metrics = [
            ("Total Tumor Volume", f"{totalTumorVol:.1f} mm3"),
            ("Total Ablation Volume", f"{totalAblationVol:.1f} mm3"),
            ("Tumor Volume Destroyed", f"{ablatedTumorVol:.1f} mm3"),
            ("Tumor Covered (Sensitivity)", f"{coverage*100:.1f} %"),
            ("Useful Ablation (Precision)", f"{precision*100:.1f} %")
        ]

        for i, (text, value) in enumerate(metrics):
            lbl = qt.QLabel(f"<b>{text}:</b> {value}")
            widget.ablationGrid.addWidget(lbl, i, 0)
            widget.ablationLabels.append(lbl)