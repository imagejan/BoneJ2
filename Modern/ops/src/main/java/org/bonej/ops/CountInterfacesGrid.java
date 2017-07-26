
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
import org.apache.commons.math3.util.FastMath;
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
		plotPlane(samplingBuilder, gridSize, sections, 0, 1);
		plotPlane(samplingBuilder, gridSize, sections, 0, 2);
		plotPlane(samplingBuilder, gridSize, sections, 1, 2);

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
		return samplingBuilder.build().map(
			s -> rotateSampler(s, unitQuaternion)).map(translateToCentre);
	}

	public static ValuePair<Vector3d, Vector3d> rotateSampler(
		final ValuePair<Vector3d, Vector3d> sampler,
		final double[] unitQuaternion)
	{
		rotateXYZ(sampler.a, unitQuaternion);
		rotateXYZ(sampler.b, unitQuaternion);
		return sampler;
	}

	public static void rotateXYZ(final Vector3d v, final double[] q) {
		final Quat4d p = new Quat4d();
		p.set(v.x, v.y, v.z, 0.0);
		final Quat4d qInv = new Quat4d();
		qInv.set(-q[0], -q[1], -q[2], q[3]);
		final Quat4d rotated = new Quat4d(q);
		rotated.mul(p);
		rotated.mul(qInv);
		v.set(rotated.x, rotated.y, rotated.z);
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
		final double gridSize, final long segments, final int dim0,
		final int dim1)
	{
		final double[] coordinates = new double[3];
		final int dim2 = orthogonalDim(dim0, dim1);
		final Vector3d normal = createNormal(dim2);
		coordinates[dim2] = -0.5 * gridSize;
		/*if (Math.random() > 0.5) {
			normal.negate();
			coordinates[dim2] = -coordinates[dim2];
		}*/
		for (double i = 0; i <= segments; i++) {
			for (double j = 0; j <= segments; j++) {
				final double offset0 = 0; // Math.random();
				final double offset1 = 0; // Math.random();
				coordinates[dim0] = (j + offset0) / segments * gridSize - 0.5 *
					gridSize;
				coordinates[dim1] = (i + offset1) / segments * gridSize - 0.5 *
					gridSize;
				samplingBuilder.add(new ValuePair<>(new Vector3d(coordinates), normal));
			}
		}
	}

	private static Vector3d createNormal(final int dim2) {
		final double[] n = new double[3];
		n[dim2] = 1.0;
		return new Vector3d(n);
	}

	private static int orthogonalDim(final int dim0, final int dim1) {
		final Set<Integer> dims = new HashSet<>(Arrays.asList(0, 1, 2));
		dims.remove(dim0);
		dims.remove(dim1);
		return dims.iterator().next();
	}

	public static boolean outOfBounds(final int[] coordinates,
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

	public static int[] toVoxelCoordinates(final Vector3d v) {
		final int[] coordinates = new int[3];
		coordinates[0] = (int) FastMath.floor(v.x);
		coordinates[1] = (int) FastMath.floor(v.y);
		coordinates[2] = (int) FastMath.floor(v.z);
		return coordinates;
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

	public static ValuePair<Double, Double> findIntersections(
		final ValuePair<Vector3d, Vector3d> sampler, final long[] bounds)
	{
		// Make min coordinates slightly negative to make vectors going parallel to
		// stack planes intersect, e.g. origin (0, 0, 0) and direction (0, 0, 1)
		final Vector3d origin = sampler.a;
		final Vector3d direction = sampler.b;
		final double minX = direction.x >= 0.0 ? 0 : bounds[0];
		final double maxX = direction.x >= 0.0 ? bounds[0] : 0;
		final double minY = direction.y >= 0.0 ? 0 : bounds[1];
		final double maxY = direction.y >= 0.0 ? bounds[1] : 0;
		final double minZ = direction.z >= 0.0 ? 0 : bounds[2];
		final double maxZ = direction.z >= 0.0 ? bounds[2] : 0;
		final double tX0 = (minX - origin.x) / direction.x;
		final double tX1 = (maxX - origin.x) / direction.x;
		final double tY0 = (minY - origin.y) / direction.y;
		final double tY1 = (maxY - origin.y) / direction.y;
		final double tZ0 = (minZ - origin.z) / direction.z;
		final double tZ1 = (maxZ - origin.z) / direction.z;
		if (tX0 > tY1 || tY0 > tX1) {
			return null;
		}
		double tMin = Math.max(tX0, tY0);
		double tMax = Math.min(tX1, tY1);
		if (tMin > tZ1 || tZ0 > tMax) {
			return null;
		}
		tMin = Math.max(tZ0, tMin);
		tMax = Math.min(tZ1, tMax);
		if (Double.isNaN(tMin) || Double.isNaN(tMax)) {
			return null;
		}

		return new ValuePair<>(tMin, tMax);
	}
}
