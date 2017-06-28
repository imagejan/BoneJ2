
package org.bonej.ops;

import static org.junit.Assert.assertEquals;

import net.imagej.ImageJ;

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
		assertEquals(7.0, samplePoint.x, 1e-15);
		assertEquals(1.0, samplePoint.y, 1e-15);
		assertEquals(-2.0, samplePoint.z, 1e-15);
	}

	@Test
	public void testFindDirection() throws Exception {
		final ValuePair<Vector3d, Vector3d> segment = new ValuePair<>(new Vector3d(1, 1, 0), new Vector3d(5, 3, 10));

		final Vector3d direction = CountInterfaces.findDirection(segment);

		assertEquals(1.0, direction.length(), 1e-15);
	}

	@Test
	public void testBox() throws Exception {

	}
}
