
package org.bonej.ops;

import static org.bonej.ops.CountInterfaces.Plane;
import static org.bonej.ops.CountInterfaces.createSamplePoints;
import static org.bonej.ops.CountInterfaces.createStackPlanes;
import static org.bonej.ops.CountInterfaces.findBounds;
import static org.bonej.ops.CountInterfaces.selectPlanes;
import static org.bonej.ops.CountInterfaces.toVoxelCoordinates;

import java.util.Arrays;
import java.util.List;
import java.util.LongSummaryStatistics;
import java.util.stream.Stream;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.ComplexType;
import net.imglib2.type.numeric.integer.LongType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesSamplingDistribution {

	private synchronized static <C extends ComplexType<C>> void sample(
		final RandomAccessibleInterval<C> interval, final Vector3d point)
	{
		final long[] coordinates = toVoxelCoordinates(point);
		final RandomAccess<C> access = interval.randomAccess(interval);
		access.setPosition(coordinates);
		final C voxel = access.get();
		voxel.setReal(voxel.getRealDouble() + 1.0);
	}

	public static void main(String... args) {
		final long width = 10;
		final long height = 10;
		final long depth = 10;
		final long voxels = width * height * depth;
		final RandomAccessibleInterval<LongType> stack = ArrayImgs.longs(width, height, depth);
		final Vector3d n = new Vector3d(1, 1, 1);
		n.normalize();
		System.out.println(n.toString());
		CountInterfaces.setSeed(0xC0FF33);
		final List<Plane> planes = createStackPlanes(stack);
		final Stream<ValuePair<Plane, Plane>> pairs = Stream.generate(
			() -> selectPlanes(planes));
		final long[] bounds = findBounds(stack);
		final Stream<Vector3d> samplePoints = createSamplePoints(pairs, 1.0,
			bounds);
		samplePoints.limit(1_000_000).sequential().forEach(point -> sample(stack, point));
		final IterableInterval<LongType> iterable = Views.flatIterable(stack);
		Cursor<LongType> cursor = iterable.cursor();
		final double[] sum = {0};
		cursor.forEachRemaining(c -> sum[0] += c.getRealDouble());
		final double totalSamples = sum[0];
		final long[] voxelSamples = new long[(int) voxels];

		cursor = iterable.cursor();
		int i = 0;
		while (cursor.hasNext()) {
			final long samples = cursor.next().get();
			voxelSamples[i] = samples;
			i++;
		}
		final double average = Arrays.stream(voxelSamples).average().orElse(0.0);
		final double stdDev = Arrays.stream(voxelSamples).mapToDouble(s -> s - average).map(s -> s * s).sum() / totalSamples;

		System.out.println("Total samples: " + totalSamples);
		System.out.println("Expected samples per voxel: " + totalSamples / voxels);
		System.out.println("Standard deviation: " + stdDev);

	}
}
