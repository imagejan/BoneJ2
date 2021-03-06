
package org.bonej.utilities;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.stream.IntStream;

import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.DefaultLinearAxis;
import net.imagej.axis.PowerAxis;
import net.imagej.units.UnitService;
import net.imglib2.IterableInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.ByteType;
import net.imglib2.type.numeric.real.DoubleType;

import org.junit.AfterClass;
import org.junit.Test;

/**
 * Unit tests for the {@link ElementUtil} class
 *
 * @author Richard Domander
 */
public class ElementUtilTest {

	public static final ImageJ IMAGE_J = new ImageJ();
	public static final UnitService unitService = IMAGE_J.context().getService(
		UnitService.class);

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}

	@Test
	public void testIsColorsBinaryFalseWhenIntervalNull() throws Exception {
		final boolean result = ElementUtil.isColorsBinary(null);

		assertFalse("A null interval should not be binary color", result);
	}

	@Test
	public void testIsColorsBinaryReturnsFalseIfIntervalEmpty() throws Exception {
		final IterableInterval<DoubleType> interval = ArrayImgs.doubles(0);

		final boolean result = ElementUtil.isColorsBinary(interval);

		assertFalse("An empty image should not be binary color", result);
	}

	@Test
	public void testIsColorsBinaryReturnsTrueForMonochrome() throws Exception {
		final IterableInterval<DoubleType> interval = ArrayImgs.doubles(2, 2);

		final boolean result = ElementUtil.isColorsBinary(interval);

		assertTrue("Monochrome image should be binary color", result);
	}

	@Test
	public void testIsColorsBinaryReturnsFalseForMulticolor() throws Exception {
		// Create a test image with many colors
		final IterableInterval<DoubleType> interval = ArrayImgs.doubles(2, 2);
		final Iterator<Integer> intIterator = IntStream.iterate(0, i -> i + 1)
			.iterator();
		interval.cursor().forEachRemaining(e -> e.setReal(intIterator.next()));

		final boolean result = ElementUtil.isColorsBinary(interval);

		assertFalse(
			"An image with more than two colours should not be binary color", result);
	}

	@Test
	public void testIsColorsBinaryBooleanTypeAssignable() throws Exception {
		final Img<BitType> img = ArrayImgs.bits(1);

		final boolean isBinary = ElementUtil.isColorsBinary(img);

		assertTrue("A BitType image should be binary", isBinary);
	}

	@Test
	public void testIsColorsBinary() throws Exception {
		// Create a test image with two colors
		final IterableInterval<DoubleType> interval = ArrayImgs.doubles(2, 2);
		final Iterator<Integer> intIterator = IntStream.iterate(0, i -> (i + 1) % 2)
			.iterator();
		interval.cursor().forEachRemaining(e -> e.setReal(intIterator.next()));

		final boolean result = ElementUtil.isColorsBinary(interval);

		assertTrue("An image with two colours should be binary color", result);
	}

	@Test
	public void testCalibratedSpatialElementSizeNullSpace() throws Exception {
		final double result = ElementUtil.calibratedSpatialElementSize(null,
			unitService);

		assertTrue("Size should be NaN when space is null", Double.isNaN(result));
	}

	@Test
	public void testCalibratedSpatialElementSizeNoSpatialAxes() throws Exception {
		final DefaultLinearAxis cAxis = new DefaultLinearAxis(Axes.CHANNEL);
		final Img<DoubleType> img = IMAGE_J.op().create().img(new int[] { 3 });
		final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", cAxis);

		final double elementSize = ElementUtil.calibratedSpatialElementSize(imgPlus,
			unitService);

		assertTrue("Element size should be NaN when there are no spatial axes",
			Double.isNaN(elementSize));
	}

	@Test
	public void testCalibratedSpatialElementSizeNonLinearAxis() throws Exception {
		final DefaultLinearAxis xAxis = new DefaultLinearAxis(Axes.X);
		final PowerAxis yAxis = new PowerAxis(Axes.Y, 2);
		final Img<DoubleType> img = IMAGE_J.op().create().img(new int[] { 10, 10 });
		final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", xAxis,
			yAxis);

		final double result = ElementUtil.calibratedSpatialElementSize(imgPlus,
			unitService);

		assertTrue("Size should be NaN when space has nonlinear axes", Double.isNaN(
			result));
	}

	@Test
	public void testCalibratedSpatialElementSizeUnitsInconvertible()
		throws Exception
	{
		final DefaultLinearAxis xAxis = new DefaultLinearAxis(Axes.X, "cm");
		final DefaultLinearAxis yAxis = new DefaultLinearAxis(Axes.Y, "");
		final Img<DoubleType> img = IMAGE_J.op().create().img(new int[] { 10, 10 });
		final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", xAxis,
			yAxis);

		final double result = ElementUtil.calibratedSpatialElementSize(imgPlus,
			unitService);

		assertEquals("Size should be 1.0 if unit inconvertible", 1.0, result, 1e-12);
	}

	@Test
	public void testCalibratedSpatialElementSizeNoUnits() throws Exception {
		final DefaultLinearAxis xAxis = new DefaultLinearAxis(Axes.X, 20.0);
		final DefaultLinearAxis yAxis = new DefaultLinearAxis(Axes.Y, 4.0);
		final DefaultLinearAxis zAxis = new DefaultLinearAxis(Axes.Z, 1.0);
		final Img<ByteType> img = ArrayImgs.bytes(1, 1, 1);
		final ImgPlus<ByteType> imgPlus = new ImgPlus<>(img, "", xAxis, yAxis,
			zAxis);

		final double elementSize = ElementUtil.calibratedSpatialElementSize(imgPlus,
			unitService);

		assertEquals("Element size is incorrect", 80.0, elementSize, 1e-12);
	}

	@Test
	public void testCalibratedSpatialElementSize() throws Exception {
		final double[][] scales = { { 20.0, 4.0, 1.0 }, { 20.0, 1.0, 4.0 }, { 4.0,
			20.0, 1.0 } };
		final String[][] units = { { "mm", "cm", "m" }, { "m", "cm", "mm" }, { "µm",
			"µm", "µm" } };
		final double[] expected = { 800_000, 0.0008, 80.0 };
		final Img<ByteType> img = ArrayImgs.bytes(1, 1, 1);
		final ImgPlus<ByteType> imgPlus = new ImgPlus<>(img);

		for (int i = 0; i < scales.length; i++) {
			final DefaultLinearAxis xAxis = new DefaultLinearAxis(Axes.X, units[i][0],
				scales[i][0]);
			final DefaultLinearAxis yAxis = new DefaultLinearAxis(Axes.Y, units[i][1],
				scales[i][1]);
			final DefaultLinearAxis zAxis = new DefaultLinearAxis(Axes.Z, units[i][2],
				scales[i][2]);
			imgPlus.setAxis(xAxis, 0);
			imgPlus.setAxis(yAxis, 1);
			imgPlus.setAxis(zAxis, 2);

			final double elementSize = ElementUtil.calibratedSpatialElementSize(
				imgPlus, unitService);

			assertEquals("Element size is incorrect", expected[i], elementSize,
				1e-12);
		}
	}
}
