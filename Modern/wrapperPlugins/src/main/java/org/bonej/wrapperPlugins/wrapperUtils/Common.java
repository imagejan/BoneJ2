
package org.bonej.wrapperPlugins.wrapperUtils;

import ij.ImagePlus;
import net.imagej.ImgPlus;
import net.imagej.axis.CalibratedAxis;
import net.imagej.ops.OpService;
import net.imglib2.img.Img;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.ComplexType;

/**
 * Miscellaneous utility methods
 *
 * @author Richard Domander
 */
public class Common {

	/** Duplicates the image without changing the title of the copy, or cropping it to the ROI */
	public static ImagePlus cleanDuplicate(final ImagePlus image) {
		image.killRoi();
		final ImagePlus copy = image.duplicate();
		image.restoreRoi();
		copy.setTitle(image.getTitle());
		return copy;
	}

	/**
	 * Converts the ImagePlus to a new ImagePlus where elements are of the given
	 * type
	 */
	public static <C extends ComplexType<C>> ImgPlus<BitType> toBitTypeImgPlus(
		OpService ops, final ImgPlus<C> imgPlus)
	{
		final Img<BitType> convertedImg = ops.convert().bit(imgPlus.getImg());
		final ImgPlus<BitType> convertedImgPlus = new ImgPlus<>(convertedImg);
		copyMetadata(imgPlus, convertedImgPlus);

		return convertedImgPlus;
	}

	/**
	 * Copies image metadata such as name, axis types and calibrations from source
	 * to target
	 */
	private static void copyMetadata(ImgPlus<?> source, ImgPlus<?> target) {
		target.setName(source.getName());

		final int dimensions = source.numDimensions();
		for (int d = 0; d < dimensions; d++) {
			final CalibratedAxis axis = source.axis(d);
			target.setAxis(axis, d);
		}
	}
}
