
package org.bonej.ops;

import static org.bonej.ops.CountInterfaces.createStackPlanes;
import static org.bonej.ops.CountInterfaces.selectPlanes;
import static org.junit.Assert.assertEquals;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.util.ValuePair;

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

		final Vector3d direction = CountInterfaces.findDirection(segment);

		assertEquals(1.0, direction.length(), 1e-15);
		assertVector(expected.x, expected.y, expected.z, direction);
	}

	@Test
	public void testSelectPlanes() throws Exception {
		final RandomAccessibleInterval<UnsignedShortType> stack = ArrayImgs
			.unsignedShorts(10, 10, 10);
		CountInterfaces.setSeed(0xC0FF33);
		final List<CountInterfaces.Plane> planes = createStackPlanes(stack);
		final long samples = 1_000_000;
		final double tolerance = 0.002;
		final double expected = samples / 6.0;

		final Collection<Long> counts = Stream.generate(() -> selectPlanes(planes)).limit(samples).collect(
				Collectors.groupingBy(Function.identity(), Collectors.counting())).values();

		//TODO test occurances of pairs with a hashing function!
		counts.forEach(c -> assertEquals(1.0, c / expected, tolerance));
	}

	private static void assertVector(final double x, final double y,
		final double z, final Vector3d actual)
	{
		assertEquals(x, actual.x, 1e-15);
		assertEquals(y, actual.y, 1e-15);
		assertEquals(z, actual.z, 1e-15);
	}
}
