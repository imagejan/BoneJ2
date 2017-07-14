
package org.bonej.ops;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.ComplexType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.util.ValuePair;

import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Quat4d;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGrid {

	private static final UnitSphereRandomVectorGenerator generator =
		new UnitSphereRandomVectorGenerator(4);
	private static final Random random = new Random(System.currentTimeMillis() +
		System.identityHashCode(new Object()));

	public static <C extends NumericType<C>> long[] findBounds(
		final RandomAccessibleInterval<C> interval)
	{
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static Stream<ValuePair<Vector3d, Vector3d>> plotSamplers(
		final double gridSize, final long sections, final long[] bounds)
	{
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> samplingBuilder = Stream
			.builder();
		final double offset0 = random.nextDouble();
		final double offset1 = random.nextDouble();
		plotPlane(samplingBuilder, bounds, gridSize, sections, 0, 1, offset0,
			offset1);
		plotPlane(samplingBuilder, bounds, gridSize, sections, 0, 2, offset0,
			offset1);
		plotPlane(samplingBuilder, bounds, gridSize, sections, 1, 2, offset0,
			offset1);

		// create a four-element vector where each element is a sampling of a normal
		// distribution. Normalize its length and you have a uniformly sampled
		// random unit quaternion which represents a uniformly sampled random
		// rotation.
		final double[] unitQuaternion = generator.nextVector();
		final Vector3d centroid = new Vector3d(bounds[0], bounds[1], bounds[2]);
		centroid.scale(0.5);
		final Function<ValuePair<Vector3d, Vector3d>, ValuePair<Vector3d, Vector3d>> translateToCentre =
			s -> {
				s.a.add(centroid);
				return s;
			};
		return samplingBuilder.build().map(s -> rotateSampler(s, unitQuaternion))
			.map(translateToCentre);
	}

	public static ValuePair<Vector3d, Vector3d> rotateSampler(
		final ValuePair<Vector3d, Vector3d> sampler,
		final double[] unitQuaternion)
	{
		return new ValuePair<>(rotateXYZ(sampler.a, unitQuaternion), rotateXYZ(
			sampler.b, unitQuaternion));
	}

	public static Vector3d rotateXYZ(Vector3d v, final double[] unitQuaternion) {
		final Quat4d q = new Quat4d(unitQuaternion);
		final Quat4d p = new Quat4d();
		p.set(v.x, v.y, v.z, 0.0);
		final Quat4d qInv = new Quat4d(unitQuaternion);
		qInv.inverse();
		final Quat4d rotated = new Quat4d();
		rotated.mul(q, p);
		rotated.mul(qInv);
		return new Vector3d(rotated.x, rotated.y, rotated.z);
	}

	private static Vector3d samplePoint(
		final ValuePair<Vector3d, Vector3d> sampler, final double t)
	{
		final Vector3d point = new Vector3d(sampler.b);
		point.scale(t);
		point.add(sampler.a);
		return point;
	}

	public static void plotPlane(
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> samplingBuilder,
		final long[] bounds, final double gridSize, final long segments,
		final int dim0, final int dim1, final double offset0,
		final double offset1)
	{
		final Set<Integer> dims = new HashSet<>(Arrays.asList(0, 1, 2));
		dims.remove(dim0);
		dims.remove(dim1);
		final int[] orthogonalDims = { dims.iterator().next() };
		final int[] planeDims = { dim0, dim1 };
		final Vector3d normal = createSparseVector(orthogonalDims, 1.0);
		final Vector3d planeShift = createSparseVector(orthogonalDims, -0.5 *
			gridSize);
		final Vector3d p0 = createSparseVector(planeDims, -0.5 * gridSize, -0.5 *
			gridSize);
		final Vector3d p1 = createSparseVector(planeDims, 0.5 * gridSize, 0.5 *
			gridSize);
		final Vector3d gridLine = new Vector3d();
		gridLine.sub(p1, p0);
		final double[] coords = new double[3];
		gridLine.get(coords);
		double prevI = 0.0;
		for (double i = 0; i <= segments; i++) {
			final double curI = i / segments * coords[dim1];
			final double iCoord = weightedAverage(offset1, curI, prevI);
			prevI = curI;
			double prevJ = 0.0;
			for (double j = 0; j <= segments; j++) {
				final double curJ = j / segments * coords[dim0];
				final double jCoord = weightedAverage(offset0, curJ, prevJ);
				prevJ = curJ;
				final Vector3d v = createSparseVector(planeDims, jCoord, iCoord);
				v.add(planeShift);
				v.add(p0);
				samplingBuilder.add(new ValuePair<>(v, normal));
			}
		}
	}

	public static double weightedAverage(final double weight, final double a,
		final double b)
	{
		return weight * a + (1.0 - weight) * b;
	}

	public static Vector3d createSparseVector(int[] dims, double... values) {
		final double[] coordinates = new double[3];
		for (int i = 0; i < dims.length; i++) {
			coordinates[dims[i]] = values[i];
		}
		return new Vector3d(coordinates);
	}

	public static boolean outOfBounds(final long[] coordinates,
		final long[] bounds)
	{
		return (coordinates[0] < 0) || (coordinates[0] >= (bounds[0])) ||
			(coordinates[1] < 0) || (coordinates[1] >= (bounds[1])) ||
			(coordinates[2] < 0) || (coordinates[2] >= (bounds[2]));
	}

	public static double findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return Math.sqrt(sumSquared);
	}

	public static long[] toVoxelCoordinates(final Vector3d v) {
		final double[] coordinates = new double[3];
		v.get(coordinates);
		final long[] voxelCoordinates = new long[3];
		for (int i = 0; i < 3; i++) {
			voxelCoordinates[i] = (long) Math.floor(coordinates[i]);
		}
		return voxelCoordinates;
	}

	public static <C extends ComplexType<C>> void sample(
		final RandomAccessibleInterval<C> interval, final long[] coordinates)
	{
		final RandomAccess<C> access = interval.randomAccess(interval);
		access.setPosition(coordinates);
		final C voxel = access.get();
		voxel.setReal(voxel.getRealDouble() + 1.0);
	}

	public static Stream<Vector3d> samplePoints(
		final Stream<ValuePair<Vector3d, Vector3d>> samplers,
		final double increment, final double gridSize)
	{
		final long iterations = (long) Math.floor(gridSize / increment) + 1;

		return samplers.flatMap(sampler -> DoubleStream.iterate(0.0, t -> t +
			increment).mapToObj(t -> CountInterfacesGrid.samplePoint(sampler, t))
			.limit(iterations));
	}
}
