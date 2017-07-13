
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
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.util.ValuePair;

import net.imglib2.view.Views;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGridTest {

	public static void main(String... args) {
		final RandomAccessibleInterval<UnsignedIntType> stack = ArrayImgs
			.unsignedInts(10, 10, 10);
		final long[] bounds = findBounds(stack);
		final double gridSize = findGridSize(bounds);
		final long rotations = 10_000_000;
		final long startTime = System.nanoTime();
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
			() -> plotSamplers(gridSize, (long) gridSize, bounds)).limit(rotations).flatMap(
				s -> s);
		final Stream<Vector3d> points = samplePoints(samplers, 1 / 2.3, gridSize);
		points.map(CountInterfacesGrid::toVoxelCoordinates).filter(
			c -> !outOfBounds(c, bounds)).forEach(c -> sample(stack, c));
		final long endTime = System.nanoTime();
		final long durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration (ms) " + durationMs);
		final IterableInterval<UnsignedIntType> iterable = Views.flatIterable(stack);
		final SummaryStatistics statistics = new SummaryStatistics();
		iterable.cursor().forEachRemaining(i -> statistics.addValue(i.getIntegerLong()));
		System.out.println(statistics.toString());
		final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		imageJ.ui().show(stack);
	}
}
