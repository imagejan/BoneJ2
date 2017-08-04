
package org.bonej.wrapperPlugins;

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
import java.util.stream.Collectors;

import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.ops.OpService;
import net.imagej.ops.special.computer.UnaryComputerOp;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

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

	@Parameter(label = "Rotations")
	private int  rotations = 1000;

	@Parameter(label = "Tolerance", required = false, persist = false)
	private double tolerance = 0.005;

	@Parameter
	private StatusService statusService;

	@Parameter
	private OpService opService;

	private UnaryFunctionOp<RandomAccessibleInterval<BitType>, Collection> milOp;

	@Override
	public void run() {
		statusService.showStatus("Anisotropy: initialising");
		final ImgPlus<BitType> bitImgPlus = Common.toBitTypeImgPlus(opService,
			inputImage);
		final List<Subspace<BitType>> subspaces = HyperstackUtils.split3DSubspaces(
			bitImgPlus).collect(Collectors.toList());
		matchOps(subspaces.get(0));
		final ExecutorService executor = Executors.newFixedThreadPool(5);
		long points = 0;
		for (Subspace<BitType> subspace : subspaces) {
			final List<Callable<Collection<?>>> tasks = new ArrayList<>(rotations);
			final long startTime = System.currentTimeMillis();
			for (int i = 0; i < rotations; i++) {
				tasks.add(() -> milOp.calculate(subspace.interval));
			}
			try {
				final List<Future<Collection<?>>> futures = executor.invokeAll(tasks);
				for (final Future<Collection<?>> future : futures) {
					final Collection<Vector3d> pointCloud = (Collection<Vector3d>) future.get();
					points += pointCloud.size();
				}
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
			} finally {
				executor.shutdown();
			}
			final long endTime = System.currentTimeMillis();
			System.out.println("Duration " + (endTime - startTime));
		}
		System.out.println("Total points: " + points);
	}

	private void matchOps(final Subspace<BitType> subspace) {
		milOp = Functions.unary(opService, MeanInterceptLengths.class, Collection.class, subspace.interval);
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
