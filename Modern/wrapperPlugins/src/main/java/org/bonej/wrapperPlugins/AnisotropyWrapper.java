
package org.bonej.wrapperPlugins;

import static org.bonej.wrapperPlugins.CommonMessages.NOT_3D_IMAGE;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_BINARY;
import static org.bonej.wrapperPlugins.CommonMessages.NO_IMAGE_OPEN;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.ops.OpService;
import net.imagej.ops.special.function.BinaryFunctionOp;
import net.imagej.ops.special.function.Functions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

import org.bonej.utilities.AxisUtils;
import org.bonej.utilities.ElementUtil;
import org.bonej.wrapperPlugins.wrapperUtils.Common;
import org.bonej.wrapperPlugins.wrapperUtils.HyperstackUtils;
import org.bonej.wrapperPlugins.wrapperUtils.HyperstackUtils.Subspace;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.ContextCommand;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * The wrapper class for calculating the degree of anisotropy using the mean
 * intercept length method.
 *
 * @author Richard Domander
 */
@Plugin(type = Command.class, menuPath = "Plugins>BoneJ>Anisotropy")
public class AnisotropyWrapper<T extends RealType<T> & NativeType<T>> extends
	ContextCommand
{

	@SuppressWarnings("unused")
	@Parameter(validater = "validateImage")
	private ImgPlus<T> inputImage;

	@Parameter(label = "Single sphere", required = false)
	private boolean singleSphere = false;

	@Parameter(label = "Vector length", initializer = "initVectorLength",
		persist = false, required = false)
	private double vectorLength;

	@Parameter(label = "Number of vectors", required = false)
	private long numVectors = 50_000;

	@Parameter(label = "Number of samples", persist = false, required = false,
		min = "1.0")
	private long samples = 2;

	@Parameter(label = "Minimum spheres", required = false, persist = false)
	private long minSpheres = 100;

	@Parameter(label = "Maximum spheres", required = false, persist = false)
	private long maxSpheres = 2_000;

	@Parameter(label = "Tolerance", required = false, persist = false)
	private double tolerance = 0.005;

	@Parameter
	private StatusService statusService;

	@Parameter
	private OpService opService;

	private BinaryFunctionOp<RandomAccessibleInterval<BitType>, Double, Long> intersectionsOp;

	@Override
	public void run() {
		statusService.showStatus("Anisotropy: initialising");
		final ImgPlus<BitType> bitImgPlus = Common.toBitTypeImgPlus(opService,
			inputImage);
		final List<Subspace<BitType>> subspaces = HyperstackUtils.split3DSubspaces(
			bitImgPlus).collect(Collectors.toList());
		matchOps(subspaces.get(0));
		subspaces.forEach(subspace -> {
			final Long intersections = intersectionsOp.calculate(subspace.interval,
				1.0);
			System.out.println(intersections);
		});
	}

	private void matchOps(final Subspace<BitType> subspace) {
		intersectionsOp = (BinaryFunctionOp) Functions.binary(opService,
			CountInterfaces.class, Long.class, subspace.interval, Double.class);
	}

	@SuppressWarnings("unused")
	private void validateImage() {
		if (inputImage == null) {
			cancel(NO_IMAGE_OPEN);
		}
		if (!ElementUtil.isColorsBinary(inputImage)) {
			cancel(NOT_BINARY);
		}
		if (AxisUtils.countSpatialDimensions(inputImage) != 3) {
			cancel(NOT_3D_IMAGE);
		}
		final int n = inputImage.numDimensions();
		final int zIndex = inputImage.dimensionIndex(Axes.Z);
		final long slices = zIndex >= 0 ? inputImage.dimension(zIndex) : 0;
		if (slices < 5) {
			cancel("Stack with at least 5 slices required");
		}
	}

	@SuppressWarnings("unused")
	private void initVectorLength() {
		vectorLength = 0.25 * getMaxDimension();
	}

	private long getMaxDimension() {
		final Optional<int[]> indices = AxisUtils.getXYZIndices(inputImage);
		if (!indices.isPresent()) {
			return 0L;
		}
		return Arrays.stream(indices.get()).mapToLong(inputImage::dimension).max()
			.orElse(0L);
	}
}
