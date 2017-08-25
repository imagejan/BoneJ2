
package org.bonej.wrapperPlugins;

import static java.util.stream.Collectors.toList;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_3D_IMAGE;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_BINARY;
import static org.bonej.wrapperPlugins.CommonMessages.NO_IMAGE_OPEN;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Stream;

import org.bonej.ops.MeanInterceptLengths;
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
import org.scijava.vecmath.Vector3d;

import net.imagej.ImgPlus;
import net.imagej.ops.OpService;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

/**
 * The wrapper class for calculating the degree of anisotropy using the mean
 * intercept length method.
 *
 * @author Richard Domander
 */
@Plugin(type = Command.class, menuPath = "Plugins>BoneJ>Anisotropy")
public class AnisotropyWrapper<T extends RealType<T> & NativeType<T>> extends ContextCommand {

	@SuppressWarnings("unused")
	@Parameter(validater = "validateImage")
	private ImgPlus<T> inputImage;

	@Parameter(label = "Rotations")
	private int rotations = 1000;

	@Parameter(label = "Tolerance", required = false, persist = false)
	private double tolerance = 0.005;

	@Parameter
	private StatusService statusService;

	@Parameter
	private OpService opService;

	private UnaryFunctionOp<RandomAccessibleInterval<BitType>, Collection<Vector3d>> milOp;

	@Override
	public void run() {
		statusService.showStatus("Anisotropy: initialising");
		final ImgPlus<BitType> bitImgPlus = Common.toBitTypeImgPlus(opService, inputImage);
		final List<Subspace<BitType>> subspaces = HyperstackUtils.split3DSubspaces(bitImgPlus).collect(toList());
		matchOps(subspaces.get(0));
		final int subspaceCount = subspaces.size();
		for (int i = 0; i < subspaceCount; i++) {
			final String ordinal = String.format("%d/%d", i + 1, subspaceCount);
			statusService.showStatus("Anisotropy: processing subspace " + ordinal);
			final Subspace<BitType> subspace = subspaces.get(i);
			final long startTime = System.currentTimeMillis();
			final Collection<Vector3d> pointCloud = subspaceMILVectors(subspace);
			if (pointCloud == null) {
				statusService.showStatus(-1, -1, "Processing failed for subspace: " + ordinal, true);
			}
			final long endTime = System.currentTimeMillis();
			System.out.println("Duration " + (endTime - startTime));
		}
	}

	private Collection<Vector3d> subspaceMILVectors(final Subspace<BitType> subspace) {
		final ExecutorService executor = Executors.newFixedThreadPool(5);
		final Collection<Vector3d> vectors = new ArrayList<>();
		try {
			final List<Future<Collection<Vector3d>>> futures = createFutures(executor, subspace);
			int calculated = 0;
			for (int i = 0; i < futures.size(); i++) {
				final Future<Collection<Vector3d>> future = futures.get(i);
				vectors.addAll(future.get());
				// TODO Fit ellipsoid to point cloud, and stop if error below
				// tolerance
				calculated++;
				statusService.showProgress(calculated, rotations);
			}
		} catch (InterruptedException | ExecutionException e) {
			return null;
		} finally {
			executor.shutdown();
		}
		return vectors;
	}

	private List<Future<Collection<Vector3d>>> createFutures(final ExecutorService executor,
			final Subspace<BitType> subspace) {
		return Stream.generate(() -> milOpCallable(subspace)).map(executor::submit).limit(rotations).collect(toList());
	}

	private Callable<Collection<Vector3d>> milOpCallable(final Subspace<BitType> subspace) {
		return () -> milOp.calculate(subspace.interval);
	}

	@SuppressWarnings("unchecked")
	private void matchOps(final Subspace<BitType> subspace) {
		milOp = (UnaryFunctionOp) Functions.unary(opService, MeanInterceptLengths.class, Collection.class,
				subspace.interval);
	}

	@SuppressWarnings("unused")
	private void validateImage() {
		if (inputImage == null) {
			cancel(NO_IMAGE_OPEN);
		}
		if (!ElementUtil.isColorsBinary(inputImage)) {
			cancel(NOT_BINARY);
		}
		// TODO Add support for 2D
		if (AxisUtils.countSpatialDimensions(inputImage) != 3) {
			cancel(NOT_3D_IMAGE);
		}
	}
}
