
package org.bonej.ops;

import static java.lang.Math.sqrt;
import static org.junit.Assert.assertEquals;

import org.junit.Test;
import org.scijava.vecmath.Vector3d;

/**
 * Tests for {@link CountInterfaces}
 *
 * @author Richard Domander
 */
public class CountInterfacesTest {

	@Test
	public void rotateUnit() throws Exception {
		final Vector3d v = new Vector3d(1.0, 0.0, 0.0);
		final Vector3d axis = new Vector3d(0.0, 0.0, 10.0);
		final Vector3d u = CountInterfaces.rotate(v, axis, -Math.PI / 2.0);
		assertVector(0.0, -1.0, 0.0, u);
	}

	@Test
	public void rotate() throws Exception {
		final Vector3d axis = new Vector3d(0.0, 1.0 / 2.0, sqrt(3.0) / 2.0);
		final Vector3d v = new Vector3d(1, -1.0, 2);
		final Vector3d u = CountInterfaces.rotate(v, axis, Math.PI / 3.0);
		assertVector((10.0 + 4.0 * sqrt(3.0)) / 8.0, (1 + 2 * sqrt(3.0)) / 8.0,
			(14 - 3 * sqrt(3.0)) / 8.0, u);

	}

	private void assertVector(final double x, final double y, final double z,
		final Vector3d v)
	{
		assertEquals(x, v.x, 1e-15);
		assertEquals(y, v.y, 1e-15);
		assertEquals(z, v.z, 1e-15);
	}

}
