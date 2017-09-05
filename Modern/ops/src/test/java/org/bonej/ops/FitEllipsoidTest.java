
package org.bonej.ops;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collection;

import net.imagej.ImageJ;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
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
public class FitEllipsoidTest {

	private static final ImageJ IMAGE_J = new ImageJ();
	private static UnaryFunctionOp<Collection<Vector3D>, Solution> fitter;
	private static Collection<Vector3D> ellipsoidPoints;

	@Test
	public void testCentre() {
		final Solution solution = fitter.calculate();
		final RealVector centre = solution.getCentre();
		assertEquals(3, centre.getDimension());
		System.out.println(centre.toString());
	}

	@BeforeClass
	public static void oneTimeSetUp() {
		// @formatter:off
		ellipsoidPoints = Arrays.asList(
				new Vector3D(1, 1, 1),
				new Vector3D(2, 1, 1),
				new Vector3D(0, 1, 1),
				new Vector3D(1, 3, 1),
				new Vector3D(1, -1, 1),
				new Vector3D(1, 1, 4),
				new Vector3D(1, 1, -2),
				new Vector3D(1, 1, 1),
				new Vector3D(1, 1, 1)
		);
		// @formatter:on
		fitter = Functions.unary(IMAGE_J.op(), FitEllipsoid.class, Solution.class,
			ellipsoidPoints);
	}

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}
}
