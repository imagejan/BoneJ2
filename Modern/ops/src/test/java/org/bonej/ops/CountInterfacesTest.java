
package org.bonej.ops;

import static org.bonej.ops.CountInterfaces.createStackPlanes;
import static org.bonej.ops.CountInterfaces.findDirection;
import static org.bonej.ops.CountInterfaces.selectPlanes;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.util.ValuePair;

import org.bonej.ops.CountInterfaces.Plane;
import org.junit.AfterClass;
import org.junit.Test;
import org.scijava.vecmath.Vector3d;

/**
 * Tests for {@link CountInterfaces}
 *
 * @author Richard Domander
 */
public class CountInterfacesTest {

	private static final ImageJ IMAGE_J = new ImageJ();

	@AfterClass
	public static void oneTimeTearDown() throws Exception {
		IMAGE_J.context().dispose();
	}

	@Test
	public void testCreateSamplePoint() throws Exception {
		final Vector3d samplePoint = CountInterfaces.createSamplePoint(new Vector3d(
			1.0, 1.0, 1.0), new Vector3d(2.0, 0.0, -1.0), 3.0);

		assertVector(7, 1, -2, samplePoint);
	}

	@Test
	public void testFindDirection() throws Exception {
		final ValuePair<Vector3d, Vector3d> segment = new ValuePair<>(new Vector3d(
			1, 1, 0), new Vector3d(5, -3, 10));
		final Vector3d expected = new Vector3d(4, -4, 10);
		expected.normalize();

		final Vector3d direction = findDirection(segment);

		assertEquals(1.0, direction.length(), 1e-15);
		assertVector(expected.x, expected.y, expected.z, direction);
	}

	@Test
	public void testPlaneDistribution() throws Exception {
		// SETUP
		final RandomAccessibleInterval<UnsignedShortType> stack = ArrayImgs
			.unsignedShorts(10, 10, 10);
		final List<Plane> planes = createStackPlanes(stack);
		final long samples = 100_000;
		final double expected = samples / numPairs(planes.size());

		// EXECUTE
		final Map<Long, Long> counts = Stream.generate(() -> selectPlanes(planes))
			.limit(samples).mapToLong(pair -> pairHash(pair, planes)).boxed().collect(
				Collectors.groupingBy(Long::longValue, Collectors.counting()));

		// VERIFY
		counts.forEach((k, v) -> assertEquals(1.0, v / expected, 0.03));
	}

	@Test
	public void testOrientationDistribution() throws Exception {
		// SETUP
		final RandomAccessibleInterval<UnsignedIntType> stack = ArrayImgs
			.unsignedInts(25, 25, 25);
		final long samples = 100_000;
		final List<Plane> planes = createStackPlanes(stack);
		final Stream<Vector3d> directions = Stream.generate(() -> selectPlanes(
			planes)).map(CountInterfaces::createSegment).map(
				CountInterfaces::findDirection);
		final double[] directionSums = new double[3];

		// EXECUTE
		directions.limit(samples).sequential().map(CountInterfaces::getOrientation)
			.forEach(o -> {
				for (int i = 0; i < 3; i++) {
					directionSums[i] += o[i];
				}
			});
		final double[] avgCoordinates = Arrays.stream(directionSums).map(c -> c /
			samples).toArray();

		// VERIFY
		final Vector3d avgDirection = new Vector3d(avgCoordinates);
		assertEquals(0.0, avgDirection.length(), 0.0075);
	}

	private long numPairs(final int items) {
		final Function<Integer, Integer> factorial = i -> IntStream.rangeClosed(1,
			i).reduce((a, b) -> a * b).orElse(0);
		return factorial.apply(items) / (2 * factorial.apply(items - 2));
	}

	private static void assertVector(final double x, final double y,
		final double z, final Vector3d actual)
	{
		assertEquals(x, actual.x, 1e-15);
		assertEquals(y, actual.y, 1e-15);
		assertEquals(z, actual.z, 1e-15);
	}

	private static long pairHash(final ValuePair<Plane, Plane> pair,
		final List<Plane> planes)
	{
		final int i = planes.indexOf(pair.a);
		final int j = planes.indexOf(pair.b);
		final int n = planes.size();
		return i <= j ? i * n + j : j * n + i;
	}
}
