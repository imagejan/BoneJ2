
package org.bonej.ops;

import static org.apache.commons.lang3.ArrayUtils.swap;
import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.outOfBounds;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;
import static org.bonej.ops.CountInterfacesGrid.sample;
import static org.bonej.ops.CountInterfacesGrid.samplePoints;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.DefaultLinearAxis;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGridTest {

	public static void plotSampling(final RandomAccessibleInterval<IntType> stack,
		final long[] bounds, final double gridSize, final long rotations)
	{
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
			() -> plotSamplers(gridSize, 100, bounds)).limit(rotations).flatMap(
				s -> s);
		final Stream<Vector3d> points = samplePoints(samplers, 1 / 2.3, gridSize);
		points.map(CountInterfacesGrid::toVoxelCoordinates).filter(
			c -> !outOfBounds(c, bounds)).forEach(c -> sample(stack, c));
	}

	/* Plot the average sampling direction */
	public static void plotSamplingOrientation(final double[][][][] orientations,
		final long[] bounds, final double gridSize, final long rotations)
	{
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
			() -> plotSamplers(gridSize, 100, bounds)).limit(rotations).flatMap(
				s -> s);

		samplers.forEach(sampler -> {
			final double increment = 1 / 2.3;
			final long iterations = (long) Math.floor(gridSize / increment) + 1;
			for (int i = 0; i < iterations; i++) {
				final double t = i * increment;
				final Vector3d point = new Vector3d(sampler.b);
				point.scale(t);
				point.add(sampler.a);
				final long[] coordinates = CountInterfacesGrid.toVoxelCoordinates(
					point);
				if (CountInterfacesGrid.outOfBounds(coordinates, bounds)) {
					continue;
				}
				orientations[(int) coordinates[0]][(int) coordinates[1]][(int) coordinates[2]][0] +=
					sampler.b.x;
				orientations[(int) coordinates[0]][(int) coordinates[1]][(int) coordinates[2]][1] +=
					sampler.b.y;
				orientations[(int) coordinates[0]][(int) coordinates[1]][(int) coordinates[2]][2] +=
					sampler.b.z;
			}
		});
	}

	public static ImgPlus<IntType> orientationsToImage(final double[][][][] orientations,
													   final long[] bounds)
	{
		final ArrayImg<IntType, IntArray> ints = ArrayImgs.ints(bounds[0],
			bounds[1], 3, bounds[2]);
		final ImgPlus<IntType> imgPlus = new ImgPlus<>(ints, "",
			new DefaultLinearAxis(Axes.X), new DefaultLinearAxis(Axes.Y),
			new DefaultLinearAxis(Axes.CHANNEL), new DefaultLinearAxis(Axes.Z));
		final RandomAccess<IntType> access = imgPlus.randomAccess();
		for (int k = 0; k < bounds[2]; k++) {
			access.setPosition(k, 3);
			for (int j = 0; j < bounds[1]; j++) {
				access.setPosition(j, 1);
				for (int i = 0; i < bounds[0]; i++) {
					access.setPosition(i, 0);
					access.setPosition(0, 2);
					final IntType element = access.get();
					element.set((int) (element.get() + Math.round(orientations[i][j][k][0])));
					access.setPosition(1, 2);
					element.set((int) (element.get() + Math.round(orientations[i][j][k][1])));
					access.setPosition(2, 2);
					element.set((int) (element.get() + Math.round(orientations[i][j][k][2])));
				}
			}
		}
		return imgPlus;
	}

	public static void printStatistics(
		final RandomAccessibleInterval<IntType> stack)
	{
		final IterableInterval<IntType> iterable = Views.flatIterable(stack);
		final SummaryStatistics statistics = new SummaryStatistics();
		iterable.cursor().forEachRemaining(i -> statistics.addValue(i.get()));
		System.out.println(statistics.toString());
	}

	public static void main(String... args) {
		// final RandomAccessibleInterval<IntType> stack = ArrayImgs.ints(100, 100,
		// 3, 100);
		final int width = 100;
		final int height = 100;
		final int depth = 100;
		final long[] bounds = { width, height, depth };
		final double gridSize = findGridSize(bounds);
		final long rotations = 100;
		final double[][][][] orientations = new double[width][height][depth][3];
		final long startTime = System.nanoTime();
		plotSamplingOrientation(orientations, bounds, gridSize, rotations);

		final long endTime = System.nanoTime();
		final long durationMs = (endTime - startTime) / 1_000_000;
		final ImgPlus<IntType> image = orientationsToImage(orientations, bounds);
		final int[] permutation = IntStream.range(0, depth).toArray();
		knuthShuffle(permutation);
		System.out.println(Arrays.toString(permutation));
		final RandomAccessibleInterval<IntType> view = Views.permuteCoordinates(image, permutation, 3);
		System.out.println("Duration (ms) " + durationMs);
		final ImageJ imageJ = new ImageJ();
		final Img<IntType> img = imageJ.op().convert().int32(Views.flatIterable(view));
		final ImgPlus<IntType> imgPlus = new ImgPlus<>(img, "",
				new DefaultLinearAxis(Axes.X), new DefaultLinearAxis(Axes.Y),
				new DefaultLinearAxis(Axes.CHANNEL), new DefaultLinearAxis(Axes.Z));
		printStatistics(imgPlus);
		imageJ.launch(args);
		imageJ.ui().show(imgPlus);
	}

	private static void knuthShuffle(final int[] indices) {
		final Random random = new Random(System.currentTimeMillis());
		final int last = indices.length - 1;
		for (int i = last; i >= 0; i--) {
			final int index = random.nextInt(i + 1);
			swap(indices, index, i);
		}
	}
}
