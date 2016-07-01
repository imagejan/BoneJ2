package org.bonej.wrapperPlugins;

import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.DefaultLinearAxis;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.DoubleType;
import org.junit.AfterClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.Matchers.anyInt;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

/**
 * Unit tests for the {@link WrapperUtils WrapperUtils} class.
 *
 * @author Richard Domander 
 */
public class WrapperUtilsTest {
    private static final ImageJ IMAGE_J = new ImageJ();

    @AfterClass
    public static void oneTimeTearDown() {
        IMAGE_J.context().dispose();
    }

    @Test
    public void testGetSizeDescription() throws Exception {
        final String[] expected = {"Size", "Area", "Volume", "Size"};
        final AxisType spatialAxis = Axes.get("Spatial", true);
        final DefaultLinearAxis axis = new DefaultLinearAxis(spatialAxis);
        final ImgPlus mockImage = mock(ImgPlus.class);
        when(mockImage.axis(anyInt())).thenReturn(axis);

        for (int i = 0; i < expected.length; i++) {
            final int dimensions = i + 1;
            when(mockImage.numDimensions()).thenReturn(dimensions);

            final String description = WrapperUtils.getSizeDescription(mockImage);

            assertTrue("Size description is incorrect", expected[i].equals(description));
        }
    }

    @Test
    public void testGetExponent() throws Exception {
        final char[] expected =
                {'\u0000', '\u00B2', '\u00B3', '\u2074', '\u2075', '\u2076', '\u2077', '\u2078', '\u2079', '\u0000'};
        final AxisType spatialAxis = Axes.get("Spatial", true);
        final DefaultLinearAxis axis = new DefaultLinearAxis(spatialAxis);
        final ImgPlus mockImage = mock(ImgPlus.class);
        when(mockImage.axis(anyInt())).thenReturn(axis);

        for (int i = 0; i < expected.length; i++) {
            final int dimensions = i + 1;
            when(mockImage.numDimensions()).thenReturn(dimensions);

            final char exponent = WrapperUtils.getExponent(mockImage);

            assertEquals("Wrong exponent character", expected[i], exponent);
        }
    }

    @Test
    public void testGetUnitHeaderReturnEmptyIfImageNull() throws Exception {
        final String result = WrapperUtils.getUnitHeader(null);

        assertTrue("Unit header should be empty", result.isEmpty());
    }

    @Test
    public void testGetUnitHeaderEmptyIfNoUnit() throws Exception {
        final DefaultLinearAxis axis = new DefaultLinearAxis(Axes.X);
        final Img<DoubleType> img = IMAGE_J.op().create().img(new int[]{10});
        final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", axis);

        final String result = WrapperUtils.getUnitHeader(imgPlus);

        assertTrue("Unit header should be empty", result.isEmpty());
    }

    @Test
    public void testGetUnitHeaderReturnEmptyIfDefaultUnitPixel() throws Exception {
        final DefaultLinearAxis axis = new DefaultLinearAxis(Axes.X, "pixel");
        final Img<DoubleType> img = IMAGE_J.op().create().img(new int[]{10});
        final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", axis);

        final String result = WrapperUtils.getUnitHeader(imgPlus);

        assertTrue("Unit header should be empty", result.isEmpty());
    }

    @Test
    public void testGetUnitHeaderReturnEmptyIfDefaultUnitUnit() throws Exception {
        final DefaultLinearAxis axis = new DefaultLinearAxis(Axes.X, "unit");
        final Img<DoubleType> img = IMAGE_J.op().create().img(new int[]{10});
        final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", axis);

        final String result = WrapperUtils.getUnitHeader(imgPlus);

        assertTrue("Unit header should be empty", result.isEmpty());
    }

    @Test
    public void testGetUnitHeader() throws Exception {
        final String unit = "mm";
        final char exponent = '³';
        final DefaultLinearAxis axis = new DefaultLinearAxis(Axes.X, unit);
        final Img<DoubleType> img = IMAGE_J.op().create().img(new int[]{10});
        final ImgPlus<DoubleType> imgPlus = new ImgPlus<>(img, "Test image", axis);

        final String result = WrapperUtils.getUnitHeader(imgPlus, exponent);

        assertEquals("Unexpected unit header", "(" + unit + exponent + ")", result);
    }
}