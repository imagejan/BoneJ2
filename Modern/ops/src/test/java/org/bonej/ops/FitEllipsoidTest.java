
package org.bonej.ops;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.DoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.bonej.ops.FitEllipsoid.Solution;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests for the {@link FitEllipsoid}
 *
 * @author Richard Domander
 */
//TODO Licence boiler plate
public class FitEllipsoidTest {

	private static final ImageJ IMAGE_J = new ImageJ();
	private static UnaryFunctionOp<Collection<Vector3D>, Solution> fitter;
	private static Solution solution;
	private static List<Vector3D> ellipsoidPoints;

	@Before
	public void setUp() {
		solution = fitter.calculate();
	}

	@Test
	public void testCentre() {
		final RealVector centre = solution.getCentre();
		// TODO void assertVector(final RealVector vector)
		assertEquals("Vector has wrong number of dimensions", 3, centre
			.getDimension());
		assertEquals("Vector has incorrect x", 1.0, centre.getEntry(0), 1e-12);
		assertEquals("Vector has incorrect y", 1.0, centre.getEntry(1), 1e-12);
		assertEquals("Vector has incorrect z", 1.0, centre.getEntry(2), 1e-12);
	}

	@Test
	public void testRadii() {
		assertEquals(1.0, solution.getA(), 1e-12);
	}

	@Test
	public void testSurface() {
		final RealMatrix surface = solution.getSurface();
		assertEquals(4, surface.getRow(0).length);
		assertTrue(MatrixUtils.isSymmetric(surface, 1e-12));
		final double a = surface.getEntry(0, 0);
		final double b = surface.getEntry(1, 1);
		final double c = surface.getEntry(2, 2);
		final double d = surface.getEntry(0, 1);
		final double e = surface.getEntry(0, 2);
		final double f = surface.getEntry(1, 2);
		final double g = surface.getEntry(0, 3);
		final double h = surface.getEntry(1, 3);
		final double i = surface.getEntry(2, 3);
		for (final Vector3D point : ellipsoidPoints) {
			final double x = point.getX();
			final double y = point.getY();
			final double z = point.getZ();
			final double polynomial =
					a * x * x +
					b * y * y +
					c * z * z +
					2 * d * x * y +
					2 * e * x * z +
					2 * f * y * z +
					2 * g * x +
					2 * h * y +
					2 * i * z;
			assertEquals(polynomial, 1.0, 1e-12);
		}
	}

	@BeforeClass
	public static void oneTimeSetUp() {
		final double a = 1;
		final double b = 2;
		final double c = 3;
		final Vector3D centre = new Vector3D(1, 1, 1);
		ellipsoidPoints = Stream.generate(
			() -> randomEllipsoidPoint(centre, a, b, c)).limit(250).collect(Collectors
				.toList());
		ellipsoidPoints.addAll(Arrays.asList(
				centre,
				new Vector3D(1 + a, 1, 1),
				new Vector3D(1 - a, 1, 1),
				new Vector3D(1, 1 + b, 1),
				new Vector3D(1, 1 - b, 1),
				new Vector3D(1, 1, 1 + c),
				new Vector3D(1, 1, 1 - c)
		));
		fitter = Functions.unary(IMAGE_J.op(), FitEllipsoid.class, Solution.class,
			ellipsoidPoints);
	}

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}

	private static Vector3D randomEllipsoidPoint(final Vector3D centre,
		final double a, final double b, final double c)
	{
		final DoubleFunction<Double> coordinate = (range) -> Math.random() * 2 *
			range - range;
		boolean inEllipsoid = false;
		double x = 0;
		double y = 0;
		double z = 0;
		while (!inEllipsoid) {
			x = coordinate.apply(a);
			y = coordinate.apply(b);
			z = coordinate.apply(c);
			inEllipsoid = ((x * x / (a * a)) + (y * y / (b * b)) + (z * z / (c *
				c))) <= 1;
		}
		return new Vector3D(x, y, z).add(centre);
	}
}
