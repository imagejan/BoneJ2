
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findBounds;
import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.outOfBounds;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;
import static org.bonej.ops.CountInterfacesGrid.sample;
import static org.bonej.ops.CountInterfacesGrid.samplePoints;

import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.type.numeric.ARGBType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.LongType;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.type.numeric.integer.UnsignedLongType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.ValuePair;

import net.imglib2.view.Views;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGridTest {

	public static void plotSampling(final RandomAccessibleInterval<IntType> stack, final long[] bounds, final double gridSize, final long rotations) {
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
				() -> plotSamplers(gridSize, 100, bounds)).limit(rotations).flatMap(
				s -> s);
		final Stream<Vector3d> points = samplePoints(samplers, 1 / 2.3, gridSize);
		points.map(CountInterfacesGrid::toVoxelCoordinates).filter(
				c -> !outOfBounds(c, bounds)).forEach(c -> sample(stack, c));
	}

	/* Plot the average sampling direction */
	public static void plotSamplingOrientation(final RandomAccessibleInterval<IntType> stack, final long[] bounds, final double gridSize, final long rotations) {
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
				() -> plotSamplers(gridSize, 100, bounds)).limit(rotations).flatMap(
				s -> s);
		final RandomAccess<IntType> access = stack.randomAccess();
		samplers.forEach(sampler -> {
			final double increment = 1 / 2.3;
			final long iterations = (long) Math.floor(gridSize / increment) + 1;
			for (int i = 0; i < iterations; i++){
				final double t = i * increment;
				final Vector3d point = new Vector3d(sampler.b);
				point.scale(t);
				point.add(sampler.a);
				final long[] coordinates = CountInterfacesGrid.toVoxelCoordinates(point);
				if (CountInterfacesGrid.outOfBounds(coordinates, bounds)) {
					continue;
				}
				access.setPosition(coordinates);
				final int value = access.get().get();
				final Vector3d n = sampler.b;
				final int r = ( (int) Math.round(n.x) & 0xff ) << 16 ;

			}
		});
	}

	public static double vectorHash(final Vector3d v) {
		return Math.round(v.z) * 5 + Math.round(v.y) * 3 + Math.round(v.x);
	}

	public static void printStatistics(final RandomAccessibleInterval<IntType> stack) {
		final IterableInterval<IntType> iterable = Views.flatIterable(stack);
		final SummaryStatistics statistics = new SummaryStatistics();
		iterable.cursor().forEachRemaining(i -> statistics.addValue(i.get()));
		System.out.println(statistics.toString());
	}

	public static void main(String... args) {
		final RandomAccessibleInterval<IntType> stack = ArrayImgs
				.ints(100, 100, 100);
		final long[] bounds = findBounds(stack);
		final double gridSize = findGridSize(bounds);
		final long rotations = 100;
		final long startTime = System.nanoTime();
		plotSamplingOrientation(stack, bounds, gridSize, rotations);
		final long endTime = System.nanoTime();
		final long durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration (ms) " + durationMs);
		printStatistics(stack);
		final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		imageJ.ui().show(stack);
	}
}
