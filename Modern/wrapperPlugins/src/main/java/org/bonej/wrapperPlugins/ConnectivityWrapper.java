
package org.bonej.wrapperPlugins;

import static org.bonej.wrapperPlugins.CommonMessages.BAD_CALIBRATION;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_3D_IMAGE;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_BINARY;
import static org.bonej.wrapperPlugins.CommonMessages.NO_IMAGE_OPEN;
import static org.scijava.ui.DialogPrompt.MessageType.INFORMATION_MESSAGE;
import static org.scijava.ui.DialogPrompt.MessageType.WARNING_MESSAGE;

import java.util.List;
import java.util.stream.Collectors;

import net.imagej.ImgPlus;
import net.imagej.ops.OpService;
import net.imagej.ops.Ops;
import net.imagej.ops.special.hybrid.Hybrids;
import net.imagej.ops.special.hybrid.UnaryHybridCF;
import net.imagej.table.DefaultColumn;
import net.imagej.table.Table;
import net.imagej.units.UnitService;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.DoubleType;

import org.bonej.utilities.AxisUtils;
import org.bonej.utilities.ElementUtil;
import org.bonej.utilities.SharedTable;
import org.bonej.wrapperPlugins.wrapperUtils.Common;
import org.bonej.wrapperPlugins.wrapperUtils.HyperstackUtils;
import org.bonej.wrapperPlugins.wrapperUtils.HyperstackUtils.Subspace;
import org.bonej.wrapperPlugins.wrapperUtils.ResultUtils;
import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.ContextCommand;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

/**
 * A wrapper UI class for the Connectivity Ops
 *
 * @author Richard Domander
 */
@Plugin(type = Command.class, menuPath = "Plugins>BoneJ>Connectivity",
	headless = true)
public class ConnectivityWrapper extends ContextCommand {

	public static final String NEGATIVE_CONNECTIVITY =
		"Connectivity is negative.\nThis usually happens if there are multiple particles or enclosed cavities.\n" +
			"Try running Purify prior to Connectivity.\n";

	@Parameter(validater = "validateImage")
	private ImgPlus<UnsignedByteType> inputImage;

	/**
	 * The connectivity results in a {@link Table}
	 * <p>
	 * Null if there are no results
	 * </p>
	 */
	@Parameter(type = ItemIO.OUTPUT, label = "BoneJ results")
	private Table<DefaultColumn<String>, String> resultsTable;

	@Parameter
	private OpService opService;

	@Parameter
	private UIService uiService;

	@Parameter
	private UnitService unitService;

	@Parameter
    private StatusService statusService;

	private UnaryHybridCF<RandomAccessibleInterval<BitType>, DoubleType> eulerCharacteristicOp;
	private UnaryHybridCF<RandomAccessibleInterval<BitType>, DoubleType> eulerCorrectionOp;

	/** A flag to avoid showing the same warning repeatedly */
	private boolean negativityWarned = false;
	/** The unit displayed in the results */
	private String unitHeader;

	@Override
	public void run() {
		statusService.showStatus("Connectivity: initialising");
		final ImgPlus<BitType> bitImgPlus = Common.toBitTypeImgPlus(opService,
			inputImage);
		final String name = inputImage.getName();
		final List<Subspace<BitType>> subspaces = HyperstackUtils.split3DSubspaces(
			bitImgPlus).collect(Collectors.toList());

		determineResultUnit();
		matchOps(subspaces.get(0).interval);
		subspaces.forEach(subspace -> {
			final String suffix = subspace.toString();
			final String label = suffix.isEmpty() ? name : name + " " + suffix;
			subspaceConnectivity(label, subspace.interval);
		});
		if (SharedTable.hasData()) {
			resultsTable = SharedTable.getTable();
		}
	}

	// region -- Helper methods --
	private void matchOps(final RandomAccessibleInterval<BitType> interval) {
        eulerCharacteristicOp = Hybrids.unaryCF(opService,
			Ops.Topology.EulerCharacteristic26NFloating.class, DoubleType.class,
                interval);
		eulerCorrectionOp = Hybrids.unaryCF(opService,
			Ops.Topology.EulerCorrection.class, DoubleType.class, interval);
	}

	private void determineResultUnit() {
		unitHeader = ResultUtils.getUnitHeader(inputImage, unitService, '³');
		if (unitHeader.isEmpty()) {
			uiService.showDialog(BAD_CALIBRATION, WARNING_MESSAGE);
		}
	}

	/** Process connectivity for one 3D subspace */
	private void subspaceConnectivity(final String label,
		final RandomAccessibleInterval<BitType> subspace)
	{
	    statusService.showStatus("Connectivity: calculating connectivity");
		final double eulerCharacteristic = eulerCharacteristicOp.calculate(subspace)
			.get();
        statusService.showStatus("Connectivity: calculating euler correction");
		final double edgeCorrection = eulerCorrectionOp.calculate(subspace).get();
		final double correctedEuler = eulerCharacteristic - edgeCorrection;
		final double connectivity = 1 - correctedEuler;
		final double connectivityDensity = calculateConnectivityDensity(subspace,
			connectivity);

		addResults(label, eulerCharacteristic, correctedEuler, connectivity,
			connectivityDensity);
	}

	private void addResults(String label, final double eulerCharacteristic,
		final double deltaEuler, final double connectivity,
		final double connectivityDensity)
	{
		if (connectivity < 0 && !negativityWarned) {
			uiService.showDialog(NEGATIVE_CONNECTIVITY, INFORMATION_MESSAGE);
			negativityWarned = true;
		}

		SharedTable.add(label, "Euler char. (χ)", eulerCharacteristic);
		SharedTable.add(label, "Corrected Euler (χ + Δχ)", deltaEuler);
		SharedTable.add(label, "Connectivity", connectivity);
		SharedTable.add(label, "Conn. density " + unitHeader, connectivityDensity);
	}

	private double calculateConnectivityDensity(
		final RandomAccessibleInterval subspace, final double connectivity)
	{
		final double elements = ((IterableInterval) subspace).size();
		final double elementSize = ElementUtil.calibratedSpatialElementSize(
			inputImage, unitService);
		return connectivity / (elements * elementSize);
	}

	@SuppressWarnings("unused")
	private void validateImage() {
		if (inputImage == null) {
			cancel(NO_IMAGE_OPEN);
			return;
		}

		if (AxisUtils.countSpatialDimensions(inputImage) != 3) {
			cancel(NOT_3D_IMAGE);
			return;
		}

		if (!ElementUtil.isColorsBinary(inputImage)) {
			cancel(NOT_BINARY);
		}
	}
	// endregion
}
