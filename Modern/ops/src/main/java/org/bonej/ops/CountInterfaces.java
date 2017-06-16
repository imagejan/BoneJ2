
package org.bonej.ops;

import net.imagej.ops.AbstractOp;
import net.imagej.ops.Op;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;

import org.scijava.ItemIO;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
@Plugin(type = Op.class)
public class CountInterfaces<B extends BooleanType<B>> extends
        AbstractOp
{
    @Parameter
    private RandomAccessibleInterval<B> interval;

    @Parameter
    private Vector3d centroid;

    @Parameter
    private Double radius;

    @Parameter(required = false)
    private Long nVectors = 50_000L;

    @Parameter(required = false, description = "Number of samples taken along the radius")
    private Long samples = 2L;

    @Parameter(type = ItemIO.OUTPUT)
    private Long interfaces;

    @Override
    public void run() {
        interfaces = 0L;
    }
}
