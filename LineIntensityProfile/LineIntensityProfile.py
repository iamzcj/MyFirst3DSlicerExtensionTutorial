import logging
import os

import slicer
import vtk
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from slicer.util import updateTableFromArray

#
# LineIntensityProfile
#

class LineIntensityProfile(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Line Intensity Profile"  # TODO: make this more human readable by adding spaces
        self.parent.categories = ["Examples"]  # TODO: set categories (folders where the module shows up in the module selector)
        self.parent.dependencies = []  # TODO: add here list of module names that this module requires
        self.parent.contributors = ["Chaojie Zheng (UIH)"]  # TODO: replace with "Firstname Lastname (Organization)"
        # TODO: update with short description of the module and a link to online module documentation
        self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#LineIntensityProfile">module documentation</a>.
"""
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

        # Additional initialization step after application startup is complete
        slicer.app.connect("startupCompleted()", registerSampleData)


#
# Register sample data sets in Sample Data module
#

def registerSampleData():
    """
    Add data sets to Sample Data module.
    """
    # It is always recommended to provide sample data for users to make it easy to try the module,
    # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

    import SampleData
    iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

    # To ensure that the source code repository remains small (can be downloaded and installed quickly)
    # it is recommended to store data sets that are larger than a few MB in a Github release.

    # LineIntensityProfile1
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category='LineIntensityProfile',
        sampleName='LineIntensityProfile1',
        # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
        # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
        thumbnailFileName=os.path.join(iconsPath, 'LineIntensityProfile1.png'),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
        fileNames='LineIntensityProfile1.nrrd',
        # Checksum to ensure file integrity. Can be computed by this command:
        #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
        checksums='SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
        # This node name will be used when the data set is loaded
        nodeNames='LineIntensityProfile1'
    )

    # LineIntensityProfile2
    SampleData.SampleDataLogic.registerCustomSampleDataSource(
        # Category and sample name displayed in Sample Data module
        category='LineIntensityProfile',
        sampleName='LineIntensityProfile2',
        thumbnailFileName=os.path.join(iconsPath, 'LineIntensityProfile2.png'),
        # Download URL and target file name
        uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
        fileNames='LineIntensityProfile2.nrrd',
        checksums='SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
        # This node name will be used when the data set is loaded
        nodeNames='LineIntensityProfile2'
    )


#
# LineIntensityProfileWidget
#

class LineIntensityProfileWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)  # needed for parameter node observation
        self.logic = None
        self._parameterNode = None
        self._updatingGUIFromParameterNode = False

    def setup(self):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath('UI/LineIntensityProfile.ui'))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = LineIntensityProfileLogic()

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
        # (in the selected parameter node).
        self.ui.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.rulerSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.rulerResolutionSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
        self.ui.tableSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.seriesSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.chartSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        # self.ui.invertOutputCheckBox.connect("toggled(bool)", self.updateParameterNodeFromGUI)

        # Buttons
        self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

    def cleanup(self):
        """
        Called when the application closes and the module widget is destroyed.
        """
        self.removeObservers()

    def enter(self):
        """
        Called each time the user opens this module.
        """
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self):
        """
        Called each time the user opens a different module.
        """
        # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
        self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    def onSceneStartClose(self, caller, event):
        """
        Called just before the scene is closed.
        """
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event):
        """
        Called just after the scene is closed.
        """
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self):
        """
        Ensure parameter node exists and observed.
        """
        # Parameter node stores all user choices in parameter values, node selections, etc.
        # so that when the scene is saved and reloaded, these settings are restored.

        self.setParameterNode(self.logic.getParameterNode())

        # Select default input nodes if nothing is selected yet to save a few clicks for the user
        if not self._parameterNode.GetNodeReference("InputVolume"):
            firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
            if firstVolumeNode:
                self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

    def setParameterNode(self, inputParameterNode):
        """
        Set and observe parameter node.
        Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
        """

        if inputParameterNode:
            self.logic.setDefaultParameters(inputParameterNode)

        # Unobserve previously selected parameter node and add an observer to the newly selected.
        # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
        # those are reflected immediately in the GUI.
        if self._parameterNode is not None:
            self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
        self._parameterNode = inputParameterNode
        if self._parameterNode is not None:
            self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

        # Initial GUI update
        self.updateGUIFromParameterNode()

    def updateGUIFromParameterNode(self, caller=None, event=None):
        """
        This method is called whenever parameter node is changed.
        The module GUI is updated to show the current state of the parameter node.
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
        self._updatingGUIFromParameterNode = True

        # Update node selectors and sliders
        self.ui.inputSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputVolume"))
        self.ui.rulerSelector.setCurrentNode(self._parameterNode.GetNodeReference("Ruler"))
        self.ui.rulerResolutionSliderWidget.value = float(self._parameterNode.GetParameter("Resolution"))
        self.ui.tableSelector.setCurrentNode(self._parameterNode.GetNodeReference("Table"))
        self.ui.seriesSelector.setCurrentNode(self._parameterNode.GetNodeReference("Series"))
        self.ui.chartSelector.setCurrentNode(self._parameterNode.GetNodeReference("Chart"))
        # self.ui.invertOutputCheckBox.checked = (self._parameterNode.GetParameter("Invert") == "true")

        # Update buttons states and tooltips
        if self._parameterNode.GetNodeReference("InputVolume") and self._parameterNode.GetNodeReference("Ruler"):
            self.ui.applyButton.toolTip = "Compute output volume"
            self.ui.applyButton.enabled = True
        else:
            self.ui.applyButton.toolTip = "Select input and output volume nodes"
            self.ui.applyButton.enabled = False

        # All the GUI updates are done
        self._updatingGUIFromParameterNode = False

    def updateParameterNodeFromGUI(self, caller=None, event=None):
        """
        This method is called when the user makes any change in the GUI.
        The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

        self._parameterNode.SetNodeReferenceID("InputVolume", self.ui.inputSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("Ruler", self.ui.rulerSelector.currentNodeID)
        self._parameterNode.SetParameter("Resolution", str(self.ui.rulerResolutionSliderWidget.value))
        self._parameterNode.SetNodeReferenceID("Table", self.ui.tableSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("Series", self.ui.seriesSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("Chart", self.ui.chartSelector.currentNodeID)
        # self._parameterNode.SetParameter("Invert", "true" if self.ui.invertOutputCheckBox.checked else "false")

        self._parameterNode.EndModify(wasModified)

    def onApplyButton(self):
        """
        Run processing when user clicks "Apply" button.
        """
        with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

            # Compute output
            self.logic.process(self.ui.inputSelector.currentNode(), self.ui.rulerSelector.currentNode(), self.ui.rulerResolutionSliderWidget.value)

            """
            # Compute inverted output (if needed)
            if self.ui.invertedOutputSelector.currentNode():
                # If additional output volume is selected then result with inverted threshold is written there
                self.logic.process(self.ui.inputSelector.currentNode(), self.ui.invertedOutputSelector.currentNode(),
                                   self.ui.rulerResolutionSliderWidget.value, not self.ui.invertOutputCheckBox.checked, showResult=False)
            """


#
# LineIntensityProfileLogic
#

class LineIntensityProfileLogic(ScriptedLoadableModuleLogic):
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

    def setDefaultParameters(self, parameterNode):
        """
        Initialize parameter node with default settings.
        """
        if not parameterNode.GetParameter("Resolution"):
            parameterNode.SetParameter("Resolution", "100.0")
        if not parameterNode.GetParameter("Invert"):
            parameterNode.SetParameter("Invert", "false")



    def process(self, inputVolume, ruler, resolution):
        """
        Run the processing algorithm.
        Can be used without GUI widget.
        :param inputVolume: volume to be thresholded
        :param outputVolume: thresholding result
        :param imageThreshold: values above/below this threshold will be set to 0
        :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
        :param showResult: show output volume in slice viewers
        """

        if not inputVolume or not ruler:
            raise ValueError("Input volume or ruler is invalid")

        import time
        startTime = time.time()
        logging.info('Processing started')

        point_ijk = [165, 60, 103]

        import numpy as np

        volumeArray = slicer.util.arrayFromVolume(inputVolume)
        value = volumeArray[point_ijk[2]][point_ijk[1]][point_ijk[0]]
        value2 = volumeArray[point_ijk[2], point_ijk[1], point_ijk[0]]
        print(f"point {point_ijk} values {value}")
        print(f"point {point_ijk} values {value2}")
        print(f"{resolution =}")

        sampleValues, distanceValues = self.probeVolume(inputVolume, ruler, resolution)

        stopTime = time.time()
        logging.info(f'Processing completed in {stopTime-startTime:.2f} seconds')

        self.plotLineIntensity(distanceValues, sampleValues)

    
    def plotLineIntensity(self, distanceValues, sampleValues):
        import numpy as np
        distanceValues = np.reshape(distanceValues, (-1, 1))
        sampleValues = np.reshape(sampleValues, (-1, 1))
        # print(f"{distanceValues =}")
        # print(f"{sampleValues =}")
        arr = np.append(distanceValues, sampleValues, axis=1)
        # print(f"{arr =}")

        parameterNode = self.getParameterNode()
        rulerNode = parameterNode.GetNodeReference("Ruler")
        rulerName = rulerNode.GetName()
        rulerColor = rulerNode.GetDisplayNode().GetSelectedColor()
        print(f"{rulerColor =}")

        # tableNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLTableNode")
        tableNode = parameterNode.GetNodeReference("Table")
        if not tableNode:
            # Save results to a new table node
            tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", f"Table_{rulerName}")
        updateTableFromArray(tableNode, arr)
        tableNode.GetTable().GetColumn(0).SetName("Distance")
        tableNode.GetTable().GetColumn(1).SetName("Intensity")
        parameterNode.SetNodeReferenceID("Table", tableNode.GetID())

        # Create plot
        # plotSeriesNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLPlotSeriesNode")
        plotSeriesNode = parameterNode.GetNodeReference("Series")
        if not plotSeriesNode:
            plotSeriesNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", f"Series_{rulerName}")
        plotSeriesNode.SetAndObserveTableNodeID(tableNode.GetID())
        plotSeriesNode.SetXColumnName("Distance")
        plotSeriesNode.SetYColumnName("Intensity")
        plotSeriesNode.SetPlotType(plotSeriesNode.PlotTypeScatter)
        # plotSeriesNode.SetColor(0, 0.6, 1.0)
        plotSeriesNode.SetColor(rulerColor)
        parameterNode.SetNodeReferenceID("Series", plotSeriesNode.GetID())

        # Create chart and add plot
        # plotChartNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLPlotChartNode")
        plotChartNode = parameterNode.GetNodeReference("Chart")
        if not plotChartNode:
            plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode", f"Chart_{rulerName}")
        plotChartNode.AddAndObservePlotSeriesNodeID(plotSeriesNode.GetID())
        # plotChartNode.YAxisRangeAutoOff()
        # plotChartNode.SetYAxisRange(0, 500000)
        parameterNode.SetNodeReferenceID("Chart", plotChartNode.GetID())

        # Show plot in layout
        slicer.modules.plots.logic().ShowChartInLayout(plotChartNode)
        """
        """


    def probeVolume(self,volumeNode,rulerNode,resolution):
        import math
        import numpy as np
        resolution = int(resolution)

        # get ruler endpoints coordinates in RAS
        p0ras = [0, 0, 0]
        rulerNode.GetNthControlPointPosition(0, p0ras)
        p0ras = [*p0ras, 1.0]
        p1ras = [0, 0, 0]
        rulerNode.GetNthControlPointPosition(1, p1ras)
        p1ras = [*p1ras, 1.0]

        print(f"{p0ras =}")
        print(f"{p1ras =}")
        length = math.sqrt((p0ras[0]-p1ras[0])**2+(p0ras[1]-p1ras[1])**2+(p0ras[2]-p1ras[2])**2)
        distanceValues = np.linspace(0, length, num=resolution+1)
        print(f"{length =}")
        # print(f"{distanceValues =}")
        
        # convert RAS to IJK coordinates of the vtkImageData
        ras2ijk = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(ras2ijk)
        p0ijk = [int(round(c)) for c in ras2ijk.MultiplyPoint(p0ras)[:3]]
        p1ijk = [int(round(c)) for c in ras2ijk.MultiplyPoint(p1ras)[:3]]
        print(f"{p0ijk =}")
        print(f"{p1ijk =}")

        # create VTK line that will be used for sampling
        line = vtk.vtkLineSource()
        line.SetResolution(resolution)
        line.SetPoint1(p0ijk[0], p0ijk[1], p0ijk[2])
        line.SetPoint2(p1ijk[0], p1ijk[1], p1ijk[2])
        
        # create VTK probe filter and sample the image
        probe = vtk.vtkProbeFilter()
        probe.SetInputConnection(line.GetOutputPort())
        probe.SetSourceData(volumeNode.GetImageData())
        probe.Update()

        sample = probe.GetOutput().GetPointData().GetArray('ImageScalars')
        nDataPoints = sample.GetNumberOfTuples()
        sampleValues = [sample.GetTuple1(i) for i in range(nDataPoints)]
        sampleValues = np.array(sampleValues)
        # for i in range(nDataPoints):
            # print(f"{i}: {sample.GetTuple1(i) =}")

        # return VTK array
        return sampleValues, distanceValues
        # return probe.GetOutput().GetPointData().GetArray('ImageScalars')


#
# LineIntensityProfileTest
#

class LineIntensityProfileTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """ Do whatever is needed to reset the state - typically a scene clear will be enough.
        """
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_LineIntensityProfile1()

    def test_LineIntensityProfile1(self):
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

        self.delayDisplay("Starting the test")

        # Get/create input data

        import SampleData
        registerSampleData()
        inputVolume = SampleData.downloadSample('LineIntensityProfile1')
        self.delayDisplay('Loaded test data set')

        inputScalarRange = inputVolume.GetImageData().GetScalarRange()
        self.assertEqual(inputScalarRange[0], 0)
        self.assertEqual(inputScalarRange[1], 695)

        outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
        threshold = 100

        # Test the module logic

        logic = LineIntensityProfileLogic()

        # Test algorithm with non-inverted threshold
        logic.process(inputVolume, outputVolume, threshold, True)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], threshold)

        # Test algorithm with inverted threshold
        logic.process(inputVolume, outputVolume, threshold, False)
        outputScalarRange = outputVolume.GetImageData().GetScalarRange()
        self.assertEqual(outputScalarRange[0], inputScalarRange[0])
        self.assertEqual(outputScalarRange[1], inputScalarRange[1])

        self.delayDisplay('Test passed')
