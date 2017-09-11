
package org.bonej.ops;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.bonej.ops.FitEllipsoid.Solution;
import org.junit.AfterClass;
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
	private static double delta = 1e-10;

	@Test
	public void testEllipsoid() {

	}

	@Test
	public void testRandomSpherePoints() {

	}

	@Test
	public void testRotatedEllipsoid() {

	}

	@Test
	public void testUnitSphere() {
		final double p = Math.cos(Math.PI / 4.0);
		final List<Vector3D> sphere = Arrays.asList(new Vector3D(1, 0, 0),
			new Vector3D(-1, 0, 0), new Vector3D(0, 1, 0), new Vector3D(0, -1, 0),
			new Vector3D(0, 0, 1), new Vector3D(0, 0, -1), new Vector3D(p, p, 0),
			new Vector3D(-p, p, 0), new Vector3D(p, -p, 0), new Vector3D(-p, -p, 0),
			new Vector3D(0, p, p), new Vector3D(0, -p, p), new Vector3D(0, p, -p),
			new Vector3D(0, -p, -p));
		final ArrayRealVector centre = new ArrayRealVector(new double[] { 0, 0,
			0 });
		// @formatter:off
		final RealVector[] eigenvectors = {
				new ArrayRealVector(new double[] { 0, 1, 0 }),
				new ArrayRealVector(new double[] { 0, 0, 1 }),
				new ArrayRealVector(new double[] { 1, 0, 0 }),
		};
		// @formatter:on
		final double[] eigenvalues = { 1, 1, 1 };

		final Solution solution = fitter.calculate(sphere);

		assertSurface(sphere, solution.getSurface());
		assertEquals(centre, solution.getCentre());
		assertArrayEquals(eigenvectors, solution.getEigenvectors());
		assertArrayEquals(eigenvalues, solution.getEigenvalues(), delta);
		assertEquals(1.0, solution.getA(), delta);
		assertEquals(1.0, solution.getB(), delta);
		assertEquals(1.0, solution.getC(), delta);
		assertTrue(solution.isEllipsoid());
	}

	public static void assertSurface(final List<Vector3D> points,
		final RealMatrix surface)
	{
		assertEquals(4, surface.getRow(0).length);
		assertTrue(MatrixUtils.isSymmetric(surface, delta));
		final double a = surface.getEntry(0, 0);
		final double b = surface.getEntry(1, 1);
		final double c = surface.getEntry(2, 2);
		final double d = surface.getEntry(0, 1);
		final double e = surface.getEntry(0, 2);
		final double f = surface.getEntry(1, 2);
		final double g = surface.getEntry(0, 3);
		final double h = surface.getEntry(1, 3);
		final double i = surface.getEntry(2, 3);
		for (final Vector3D point : points) {
			final double x = point.getX();
			final double y = point.getY();
			final double z = point.getZ();
			final double polynomial = a * x * x + b * y * y + c * z * z + 2 * d * x *
				y + 2 * e * x * z + 2 * f * y * z + 2 * g * x + 2 * h * y + 2 * i * z;
			assertEquals(1.0, polynomial, delta);
		}
	}

	@BeforeClass
	public static void oneTimeSetUp() {
		final Collection<Vector3D> points = Stream.generate(() -> new Vector3D(0, 0,
			0)).limit(FitEllipsoid.SURFACE_TERMS).collect(Collectors.toList());
		fitter = Functions.unary(IMAGE_J.op(), FitEllipsoid.class, Solution.class,
			points);
	}

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}
}
