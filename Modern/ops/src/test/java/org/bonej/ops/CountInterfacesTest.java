
package org.bonej.ops;

import net.imagej.ImageJ;

import org.junit.AfterClass;

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

}
