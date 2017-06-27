
package org.bonej.wrapperPlugins;

import net.imagej.ImageJ;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.LongArray;
import net.imglib2.type.logic.BitType;
import org.bonej.ops.CountInterfaces;

/**
 * A main class for quickly testing the wrapper plugins
 *
 * @author Richard Domander
 */
public class BoneJMain {

	public static void main(String... args) {
		final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		final ArrayImg<BitType, LongArray> bits = ArrayImgs.bits(11, 11, 11);
		imageJ.op().run(CountInterfaces.class, bits, 0.5);
	}
}
