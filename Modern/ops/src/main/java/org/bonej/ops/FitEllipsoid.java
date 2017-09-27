
package org.bonej.ops;

import java.util.Arrays;
import java.util.Collection;
import java.util.Optional;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.stream.IntStream;

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.scijava.plugin.Plugin;

/**
 * An op that tries to fit an ellipsoid on a set of points.
 * <p>
 * The op first solves the quadric surface that fits the points by minimising
 * the distance by least squares fitting. This surface is found by solving a
 * nine term polynomial. It may or may not be an ellipsoid. After the solution's
 * been found, the op calculates the properties of the ellipsoid such as center
 * point, radii and eigenvectors that describe the orientation.
 * </p>
 * <p>
 * The op is based on the <a href=
 * "https://www.researchgate.net/publication/4070857_Least_squares_ellipsoid_specific_fitting">algorithm</a>
 * of Li &amp; Griffiths (DOI: 10.1109/GMAP.2004.1290055), and the
 * implementations of Yury Petrov &amp; "KalebKE".
 * </p>
 * 
 * @author Richard Domander
 */
@Plugin(type = Op.class)
public class FitEllipsoid extends
	AbstractUnaryFunctionOp<Collection<Vector3D>, Optional<FitEllipsoid.Solution>>
	implements Contingent
{

	/**
	 * Number of terms in the surface polynomial that needs to be solved.
	 * <p>
	 * Due the number of terms, we also need at least 9 points in the input to
	 * solve the equation.
	 * </p>
	 * 
	 * @see #solveSurface(Collection)
	 */
	public static final int SURFACE_TERMS = 9;

	/**
	 * Minimum value for an eigenvalue to be considered non-zero.
	 * 
	 * @see #isEllipsoid(double[])
	 */
	private static final double EIGENVALUE_TOLERANCE = 1e-10;

	/**
	 * Tries to solve an ellipsoid that fits the points.
	 *
	 * @param points a point cloud in 3D space.
	 * @return an {@link Optional} containing the properties of the solved
	 *         ellipsoid, or empty if an ellipsoidal solution was not found.
	 */
	@Override
	public Optional<Solution> calculate(final Collection<Vector3D> points) {
		final RealVector polynomial = solveSurface(points);
		final RealMatrix surface = toAlgebraicMatrix(polynomial);
		final RealVector center = solveCenter(surface);
		final RealMatrix translatedSurface = translateToCenter(surface, center);
		final EigenDecomposition decomposition = solveEigenvectors(
			translatedSurface);
		final double[] eigenvalues = decomposition.getRealEigenvalues();
		if (!isEllipsoid(eigenvalues)) {
			return Optional.empty();
		}
		final double[] radii = calculateRadii(eigenvalues);
		return Optional.of(new Solution(surface, center, decomposition, radii));
	}

	/** {@inheritDoc} */
	@Override
	public boolean conforms() {
		return in().size() >= SURFACE_TERMS;
	}

	// region -- Helper methods --

	/**
	 * Converts eigenvalues to the radii of an ellipsoid.
	 *
	 * @param eigenvalues positive eigenvalues of a quadric.
	 * @return the radii in ascending order.
	 */
	private double[] calculateRadii(final double[] eigenvalues) {
		final DoubleUnaryOperator toRadius = x -> FastMath.sqrt(1.0 / x);
		return Arrays.stream(eigenvalues).map(toRadius).sorted().toArray();
	}

	/**
	 * Creates a design matrix out of a collection of points.
	 * <p>
	 * The design matrix is used in the least squares fitting of a quadric
	 * surface.
	 * </p>
	 *
	 * @see #solveSurface(Collection)
	 * @param pointCloud points in a 3D space (size >= {@link #SURFACE_TERMS}).
	 * @return a [n][9] matrix of real values.
	 */
	private static RealMatrix createDesignMatrix(
		final Collection<Vector3D> pointCloud)
	{
		final Function<Vector3D, double[]> toRow = (point) -> {
			final double x = point.getX();
			final double y = point.getY();
			final double z = point.getZ();
			return new double[] { x * x, y * y, z * z, 2 * x * y, 2 * x * z, 2 * y *
				z, 2 * x, 2 * y, 2 * z };
		};
		final double[][] data = pointCloud.stream().map(toRow).toArray(
			double[][]::new);
		return new Array2DRowRealMatrix(data);
	}

	/**
	 * Determines if the eigenvalues belong to an ellipsoid.
	 * <p>
	 * The signs of the eigenvalues determine the type of a quadric. If they are
	 * all positive, it is an ellipsoid; two positive and one negative gives a
	 * hyperboloid of one sheet; one positive and two negative gives a hyperboloid
	 * of two sheets. If one or more eigenvalues vanish, we have a degenerate case
	 * such as a paraboloid,or a cylinder or even a pair of planes.
	 * </p>
	 *
	 * @param eigenvalues eigenvalues of a quadric surface
	 * @return true if the surface is an ellipsoid, false otherwise.
	 */
	private static boolean isEllipsoid(final double[] eigenvalues) {
		return Arrays.stream(eigenvalues).allMatch(x -> x > EIGENVALUE_TOLERANCE);
	}

	/**
	 * Finds the center point of a quadric surface.
	 *
	 * @param surface the general equation of a surface in algebraic matrix form.
	 * @see #toAlgebraicMatrix(RealVector)
	 * @return the 3D center point in a vector.
	 */
	private static RealVector solveCenter(final RealMatrix surface) {
		final RealMatrix subMatrix = surface.getSubMatrix(0, 2, 0, 2)
			.scalarMultiply(-1.0);
		final RealMatrix subInverse = new SingularValueDecomposition(subMatrix)
			.getSolver().getInverse();
		// The {x,y,z} translation (from origin) part of the matrix
		final RealVector translation = surface.getRowVector(3).getSubVector(0, 3);
		return subInverse.operate(translation);
	}

	/**
	 * Creates the {@link EigenDecomposition} of the given quadric surface.
	 *
	 * @param surface the general equation of a surface in algebraic matrix form.
	 * @see #toAlgebraicMatrix(RealVector)
	 * @return a new copy of the eigen decomposition.
	 */
	private static EigenDecomposition solveEigenvectors(
		final RealMatrix surface)
	{
		final double scalar = -1.0 / surface.getEntry(3, 3);
		final RealMatrix eigenMatrix = surface.getSubMatrix(0, 2, 0, 2)
			.scalarMultiply(scalar);
		return new EigenDecomposition(eigenMatrix);
	}

	/**
	 * Solves the quadratic surface that best fits the given points.
	 * <p>
	 * The equation solved is the polynomial Ax<sup>2</sup> + By<sup>2</sup> +
	 * Cz<sup>2</sup> + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz, i.e. the general
	 * equation of a quadric surface.
	 * </p>
	 *
	 * @param pointCloud A collection of points in a 3D space.
	 * @return the solution vector of the surface. The solution may not form an
	 *         ellipsoid.
	 */
	private static RealVector solveSurface(
		final Collection<Vector3D> pointCloud)
	{
		final RealMatrix d = createDesignMatrix(pointCloud);
		final RealMatrix dT = d.transpose();
		final RealMatrix dTDInverse = new SingularValueDecomposition(dT.multiply(d))
			.getSolver().getInverse();
		final ArrayRealVector ones = new ArrayRealVector(pointCloud.size(), 1);
		return dTDInverse.operate(dT.operate(ones));
	}

	/**
	 * Creates a matrix out of a quadric surface solution vector.
	 *
	 * @see #solveSurface(Collection)
	 * @return a matrix representing the polynomial solution vector in an
	 *         algebraic form.
	 */
	private static RealMatrix toAlgebraicMatrix(final RealVector v) {
		// I'm not a clever man, so I'm using named placeholder variables to better
		// follow the array assignment
		final double a = v.getEntry(0);
		final double b = v.getEntry(1);
		final double c = v.getEntry(2);
		final double d = v.getEntry(3);
		final double e = v.getEntry(4);
		final double f = v.getEntry(5);
		final double g = v.getEntry(6);
		final double h = v.getEntry(7);
		final double i = v.getEntry(8);
		return new Array2DRowRealMatrix(new double[][] { { a, d, e, g }, { d, b, f,
			h }, { e, f, c, i }, { g, h, i, -1 } });
	}

	/**
	 * Translates the quadratic surface equation to the center point.
	 *
	 * @see #toAlgebraicMatrix(RealVector)
	 * @param surface the general equation of a surface in algebraic matrix form.
	 * @param center a 3D center point.
	 * @return translated surface matrix.
	 */
	private static RealMatrix translateToCenter(final RealMatrix surface,
		final RealVector center)
	{
		final double x = center.getEntry(0);
		final double y = center.getEntry(1);
		final double z = center.getEntry(2);
		final Array2DRowRealMatrix t = new Array2DRowRealMatrix(new double[][] { {
			1, 0, 0, x }, { 0, 1, 0, y }, { 0, 0, 1, z }, { 0, 0, 0, 1 } });
		return t.transpose().multiply(surface).multiply(t);
	}
	// endregion

	/** Stores the properties of an ellipsoid solved by the op. */
	// region -- Helper classes --
	public static class Solution {

		private final RealMatrix surface;
		private final RealVector center;
		private final EigenDecomposition decomposition;
		private final double[] radii;

		private Solution(final RealMatrix surface, final RealVector center,
			final EigenDecomposition decomposition, final double[] radii)
		{
			this.surface = surface;
			this.center = center;
			this.decomposition = decomposition;
			this.radii = radii;
		}

		/**
		 * Gets the radius of the first axis of the ellipsoid.
		 *
		 * @return a real radius of the ellipsoid.
		 */
		public double getA() {
			return radii[0];
		}

		/**
		 * Gets the radius of the second axis of the ellipsoid.
		 *
		 * @return a real radius of the ellipsoid.
		 */
		public double getB() {
			return radii[1];
		}

		/**
		 * Gets the radius of the third axis of the ellipsoid.
		 *
		 * @return a real radius of the ellipsoid.
		 */
		public double getC() {
			return radii[2];
		}

		/**
		 * Gets a copy of the center of the ellipsoid.
		 *
		 * @return A copy of the center vector.
		 */
		public RealVector getCenter() {
			return center.copy();
		}

		/**
		 * Gets a copy of eigenvalues of the ellipsoid surface.
		 *
		 * @return copy of the real parts of the eigenvalues.
		 */
		public double[] getEigenvalues() {
			return decomposition.getRealEigenvalues();
		}

		/**
		 * Get a copy of the eigenvectors of the ellipsoid surface.
		 *
		 * @return eigenvector array.
		 */
		public RealVector[] getEigenvectors() {
			return new RealVector[] { decomposition.getEigenvector(0), decomposition
				.getEigenvector(1), decomposition.getEigenvector(2) };
		}

		/**
		 * Gets a copy of the radii of the ellipsoid.
		 *
		 * @return an array of radii {a, b, c}.
		 */
		public double[] getRadii() {
			return Arrays.copyOf(radii, radii.length);
		}

		/**
		 * Gets the semiaxes of the ellipsoid.
		 *
		 * @return an array of vectors {u, v, w}.
		 */
		public RealVector[] getSemiaxes() {
			final RealVector[] vectors = getEigenvectors();
			return IntStream.range(0, 3).mapToObj(i -> vectors[i].mapMultiply(
				radii[i])).toArray(RealVector[]::new);
		}

		/**
		 * Gets the algebraic form of the ellipsoid surface.
		 *
		 * @return matrix of the surface terms.
		 */
		public RealMatrix getSurface() {
			return surface.copy();
		}

	}
	// endregion
}
