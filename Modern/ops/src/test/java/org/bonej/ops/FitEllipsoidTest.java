
package org.bonej.ops;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Stream.generate;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.function.DoubleBinaryOperator;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.bonej.ops.FitEllipsoid.Solution;
import org.hamcrest.CoreMatchers;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;
import org.scijava.vecmath.Quat4d;

/**
 * Tests for the {@link FitEllipsoid}
 * <p>
 * What's missing is a test for random points within an ellipsoid. There the
 * fitting is too unpredictable to check for expected values with any
 * meaningful margin of error. The solved surface may not even be an ellipsoid.
 * There may be an additional constraint that ensures that the generated points
 * have an ellipsoidal solution, but I don't know of such.
 * </p>
 * 
 * @author Richard Domander
 */
public class FitEllipsoidTest {

	private static final ImageJ IMAGE_J = new ImageJ();
	private static final double DELTA = 1e-10;
	private static final long SEED = 0xC0FF33;
	private static final UnitSphereRandomVectorGenerator quaternionGenerator =
		new UnitSphereRandomVectorGenerator(4, new MersenneTwister(SEED));
	private static final Random random = new Random(SEED);
	private static UnaryFunctionOp<Collection<Vector3D>, Optional<Solution>> fitter;
	@Rule
	public final ExpectedException exception = ExpectedException.none();

	/**
	 * Tests ellipsoid fitting on an ellipsoidal band of points around certain
	 * centre point
	 */
	@Test
	public void testEllipsoidBand() {
		// SETUP
		final double a = 1.0;
		final double b = 2.0;
		final double c = 3.0;
		final Vector3D centre = new Vector3D(1, 1, 1);
		final List<Vector3D> points = generateEllipsoidPoints(1_000, a, b, c, 0.05,
			centre);

		// EXECUTE
		final Optional<Solution> result = fitter.calculate(points);

		// VERIFY
		assertTrue("Fitting should have found an ellipsoid", result.isPresent());
		final Solution solution = result.get();
		assertArrayEquals("Ellipsoid centre point is not within tolerance", centre
			.toArray(), solution.getCenter().toArray(), 0.05);
		assertEquals("Ellipsoid radius is not within tolerance", a, solution.getA(),
			0.025);
		assertEquals("Ellipsoid radius is not within tolerance", b, solution.getB(),
			0.025);
		assertEquals("Ellipsoid radius is not within tolerance", c, solution.getC(),
			0.025);
	}

	@Test
	public void testKnownUnitSphere() {
		final double[] centreCoords = { 0, 0, 0 };
		final double p = Math.cos(Math.PI / 4.0);
		final List<Vector3D> sphere = Stream.of(new Vector3D(1, 0, 0), new Vector3D(
			-1, 0, 0), new Vector3D(0, 1, 0), new Vector3D(0, -1, 0), new Vector3D(0,
				0, 1), new Vector3D(0, 0, -1), new Vector3D(p, p, 0), new Vector3D(-p,
					p, 0), new Vector3D(p, -p, 0), new Vector3D(-p, -p, 0), new Vector3D(
						0, p, p), new Vector3D(0, -p, p), new Vector3D(0, p, -p),
			new Vector3D(0, -p, -p)).collect(toList());
		// @formatter:off
		final RealVector[] eigenvectors = {
				new ArrayRealVector(new double[] { 0, 1, 0 }),
				new ArrayRealVector(new double[] { 0, 0, 1 }),
				new ArrayRealVector(new double[] { 1, 0, 0 }),
		};
		// @formatter:on
		final double[] eigenvalues = { 1, 1, 1 };

		final Optional<Solution> result = fitter.calculate(sphere);

		assertTrue("Fitting should have found an ellipsoid", result.isPresent());
		final Solution solution = result.get();
		assertSurface(sphere, solution.getSurface());
		assertArrayEquals("Centre point is incorrect", centreCoords, solution
			.getCenter().toArray(), DELTA);
		assertArrayEquals("Eigenvectors are incorrect", eigenvectors, solution
			.getEigenvectors());
		assertArrayEquals("Eigenvalues are incorrect", eigenvalues, solution
			.getEigenvalues(), DELTA);
		assertEquals("Radius A is incorrect", 1.0, solution.getA(), DELTA);
		assertEquals("Radius B is incorrect", 1.0, solution.getB(), DELTA);
		assertEquals("Radius C is incorrect", 1.0, solution.getC(), DELTA);
		assertArrayEquals(new double[]{1, 1, 1}, solution.getRadii(), DELTA);
		// Since radii = 1, eigenvectors are equal to the semiaxes
		final RealVector[] semiaxes = solution.getSemiaxes();
		for (int i = 0; i < 3; i++) {
			assertArrayEquals(eigenvectors[i].toArray(), semiaxes[i].toArray(), DELTA);
		}
	}

	@Test
	public void testPlanePair() {
		final List<Vector3D> planePair = Stream.generate(() -> {
			final double x = random.nextDouble() * 10;
			final double y = random.nextDouble() * 10;
			final double z = random.nextDouble() > 0.5 ? 100 : 0;
			return new Vector3D(x, y, z);
		}).limit(100).collect(toList());

		final Optional<Solution> result = fitter.calculate(planePair);

		assertFalse("Fitting should not have succeeded", result.isPresent());
	}

	@Test
	public void testRotatedEllipsoid() {
		// SETUP
		final double a = 1.0;
		final double b = 2.0;
		final double c = 3.0;
		final long n = 1_000;
		final Vector3D centre = new Vector3D(0, 0, 0);
		final double[] q = quaternionGenerator.nextVector();
		final List<Vector3D> points = generateEllipsoidPoints(n, a, b, c, 0.0,
			centre);
		assertTrue("Sanity check failed: point(s) not on ellipsoid surface", points
			.stream().allMatch(p -> isEllipsoidPoint(p, a, b, c)));
		final List<Vector3D> rotated = points.stream().map(p -> rotate(p, q))
			.collect(toList());
		for (int i = 0; i < n; i++) {
			assertEquals("Sanity check failed: rotation incorrect", points.get(i)
				.getNormSq(), rotated.get(i).getNormSq(), DELTA);
		}

		// EXECUTE
		final Optional<Solution> result = fitter.calculate(points);

		// VERIFY
		assertTrue("Fitting should have found an ellipsoid", result.isPresent());
		final Solution solution = result.get();
		assertArrayEquals("Ellipsoid centre point is not within tolerance", centre
			.toArray(), solution.getCenter().toArray(), 0.05);
		assertEquals("Ellipsoid radius is not within tolerance", a, solution.getA(),
			0.025);
		assertEquals("Ellipsoid radius is not within tolerance", b, solution.getB(),
			0.025);
		assertEquals("Ellipsoid radius is not within tolerance", c, solution.getC(),
			0.025);
		assertArrayEquals("Radii are not in order", new double[]{a, b, c}, solution.getRadii(), 0.025);
	}

	@Test
	public void testTooFewPoints() {
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage(CoreMatchers.containsString(
			"Inputs do not conform to op rules"));
		final List<Vector3D> fewPoints = generate(() -> new Vector3D(0, 0, 0))
			.limit(FitEllipsoid.SURFACE_TERMS - 1).collect(toList());
		IMAGE_J.op().run(FitEllipsoid.class, fewPoints);
	}

	@SuppressWarnings("unchecked")
	@BeforeClass
	public static void oneTimeSetUp() {
		final Collection<Vector3D> mockPoints = generate(() -> new Vector3D(0, 0,
			0)).limit(FitEllipsoid.SURFACE_TERMS).collect(toList());
		fitter = (UnaryFunctionOp) Functions.unary(IMAGE_J.op(), FitEllipsoid.class,
			Optional.class, mockPoints);
	}

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}

	private static void assertSurface(final List<Vector3D> points,
		final RealMatrix surface)
	{
		assertEquals(4, surface.getRow(0).length);
		assertTrue(MatrixUtils.isSymmetric(surface, DELTA));
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
			assertEquals("Polynomial does not solve surface equation", 1.0,
				polynomial, DELTA);
		}
	}

	/**
	 * Generates random uniformly distributed points on an ellipsoid surface.
	 * <p>
	 * NB Radii should be 0 &lt; a &lt; b &lt; c.
	 * </p>
	 * <p>
	 * NB Distribution is not uniform if variance &neq;
	 * </p>
	 * <p>
	 * In essence the method first generates uniformly distributed points on a
	 * sphere, which are randomly discarded with probability p(x, y, z). The
	 * probability is the smaller the larger the surface area of the ellipsoid
	 * would be in relation to the sphere around the point. This counteracts the
	 * tendency of points to cluster around the poles of the shortest axis.
	 * Finally the points are mapped to the ellipsoid surface.
	 * </p>
	 * <p>
	 * <a href=
	 * "https://math.stackexchange.com/questions/973101/how-to-generate-points-uniformly-distributed-on-the-surface-of-an-ellipsoid">Implementation
	 * explained in detail.</a>
	 * </p>
	 *
	 * @param n number of points generated.
	 * @param a smallest radius of ellipsoid
	 * @param b second radius of ellipsoid
	 * @param c largest radius of ellipsoid
	 * @param variance percentage radii can vary to create a band around the
	 *          surface
	 * @param centre centre point of the ellipsoid
	 */
	private static List<Vector3D> generateEllipsoidPoints(final long n,
		final double a, final double b, final double c, final double variance,
		final Vector3D centre)
	{
		final Vector3D u = new Vector3D(1, 0, 0);
		final double muMax = b * c;
		// Probability function to keep a sphere point
		final Predicate<Vector3D> pGt = v -> random.nextDouble() <= muFactor(v, a,
			b, c) / muMax;
		// Mapping function from sphere to ellipsoid
		final Function<Vector3D, Vector3D> toEllipsoid = v -> new Vector3D(a * v
			.getX(), b * v.getY(), c * v.getZ());
		// Random scaling to put the point on a band around the ellipsoid surface
		final Function<Vector3D, Vector3D> scaling = v -> {
			final double scalar = (2 * random.nextDouble() - 1) * variance + 1.0;
			return v.scalarMultiply(scalar);
		};
		return generate(() -> randomRotation(u)).filter(pGt).limit(n).map(
			toEllipsoid).map(scaling).map(v -> v.add(centre)).collect(toList());
	}

	private static boolean isEllipsoidPoint(final Vector3D p, final double a,
		final double b, final double c)
	{
		final DoubleBinaryOperator f = (x, y) -> (x * x) / (y * y);
		final double solution = f.applyAsDouble(p.getX(), a) + f.applyAsDouble(p
			.getY(), b) + f.applyAsDouble(p.getZ(), c);
		return Math.abs(1.0 - solution) <= DELTA;
	}

	private static double muFactor(final Vector3D point, final double a,
		final double b, final double c)
	{
		final double u = a * c * point.getY();
		final double v = a * b * point.getZ();
		final double w = b * c * point.getX();
		return Math.sqrt(u * u + v * v + w * w);
	}

	private static Vector3D randomRotation(final Vector3D v) {
		final double[] quaternion = quaternionGenerator.nextVector();
		return rotate(v, quaternion);
	}

	private static Vector3D rotate(final Vector3D v, final double[] quaternion) {
		final Quat4d rotated = new Quat4d(quaternion);
		final Quat4d p = new Quat4d();
		p.set(v.getX(), v.getY(), v.getZ(), 0);
		rotated.mul(p);
		final Quat4d qInv = new Quat4d(quaternion);
		qInv.inverse();
		rotated.mul(qInv);
		return new Vector3D(rotated.x, rotated.y, rotated.z);
	}
}
