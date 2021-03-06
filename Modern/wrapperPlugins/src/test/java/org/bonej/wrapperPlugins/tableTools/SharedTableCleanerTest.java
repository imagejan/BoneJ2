
package org.bonej.wrapperPlugins.tableTools;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import net.imagej.ImageJ;

import org.bonej.utilities.SharedTable;
import org.junit.AfterClass;
import org.junit.Test;

/**
 * Tests for {@link SharedTableCleaner}
 *
 * @author Richard Domander
 */
public class SharedTableCleanerTest {

	private static final ImageJ IMAGE_J = new ImageJ();

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}

	@Test
	public void testRun() throws Exception {
		// SETUP
		SharedTable.add("Label", "Header", "Value");
		assertTrue("Sanity check failed, no data in SharedTable", SharedTable
			.hasData());

		// EXECUTE
		IMAGE_J.command().run(SharedTableCleaner.class, true).get();

		// VERIFY
		assertFalse("Table should have no data", SharedTable.hasData());
	}

}
