
package org.bonej.wrapperPlugins.wrapperUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
import java.util.stream.Stream;
import java.util.stream.Stream.Builder;

import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.TypedAxis;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.LongType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

import org.bonej.wrapperPlugins.wrapperUtils.HyperstackUtils.Subspace.HyperAxisMeta;

/**
 * A static class containing utilities for splitting an n-dimensional hyperstack
 * into arbitrary subspaces
 * <p>
 * Doesn't copy the subspaces, rather provides {@link Subspace} objects, which
 * contain a {@link RandomAccessibleInterval} that can be used to traverse a
 * certain subspace. Each {@link Subspace} also contains metadata, which locates
 * the subspace in the hyperspace.
 * </p>
 *
 * @author Richard Domander
 * @implNote This code is hacky spaghetti, and it's design is all over the
 *           place. It works for now, but a refactor is in place if/when
 *           {@link ImgPlus} API changes or it's replaced by another metadata
 *           rich class
 * @implNote If you want to split a hyperspace into {X, Y, T} subspaces, and the
 *           hyperspace has more than one time dimension, *all* of the time
 *           dimensions will be lumped into {X, Y, T1, T2, .. TN} subspaces.
 *           That is, instead of {X, Y, T1}, {X, Y, T2}, .. {X, Y, TN}.
 */
public class HyperstackUtils {

	private HyperstackUtils() {}

	/**
	 * Splits the hyperstack into {X, Y, Z} subspaces
	 * 
	 * @see #splitSubspaces(ImgPlus, List)
	 */
	public static <T extends RealType<T> & NativeType<T>> Stream<Subspace<T>>
		split3DSubspaces(final ImgPlus<T> hyperStack)
	{
		return splitSubspaces(hyperStack, Arrays.asList(Axes.X, Axes.Y, Axes.Z));
	}

	/**
	 * Splits the hyperstack into subspaces defined by the given axes
	 * <p>
	 * If all the given axis types are not found in the hyperstack, gives
	 * subspaces of the found types. If none of the types are found, returns an
	 * empty stream. For example, if you want to split a {X, Y, C, T} hyperstack
	 * into {X, Y, Z}, returns all the {X, Y} subspaces.
	 * </p>
	 * <p>
	 * NB Assumes that the given {@link ImgPlus} has the necessary metadata, i.e.
	 * its {@link CalibratedAxis} have {@link AxisType}.
	 * </p>
	 *
	 * @param hyperStack An n-dimensional image
	 * @param subspaceTypes The types of the axis in the desired subspace
	 * @return A stream of all the subspaces found
	 */
	public static <T extends RealType<T> & NativeType<T>> Stream<Subspace<T>>
		splitSubspaces(final ImgPlus<T> hyperStack, List<AxisType> subspaceTypes)
	{
		if (subspaceTypes == null || subspaceTypes.isEmpty()) {
			return Stream.empty();
		}
		final Builder<Subspace<T>> builder = Stream.builder();
		final int[] splitIndices = findSplitAxisIndices(hyperStack, subspaceTypes);
		final long[] typeSubscripts = mapTypeSubscripts(hyperStack, splitIndices);
		final int numSplits = splitIndices.length;
		final HyperAxisMeta[] subspaceMeta = new HyperAxisMeta[numSplits];
		final List<ValuePair<IntType, LongType>> split = new ArrayList<>();
		splitDims(hyperStack, splitIndices, numSplits - 1, subspaceMeta, split,
			typeSubscripts, builder);
		return builder.build();
	}

	// region -- Helper methods --

	/**
	 * Maps the type subscripts of axes in the hyperstack
	 * <p>
	 * If the hyperstack has multiple axes of the same type, e.g. more than one
	 * time-axis, the subscripts are used to tell them apart. The default
	 * subscript for each axis type is one.
	 * </p>
	 * 
	 * @param hyperStack the space to be split
	 * @param splitIndices Indices of the axes
	 * @return An array of subscripts where [i] is the subscript for axis at
	 *         splitIndices[i]
	 */
	private static <T extends RealType<T> & NativeType<T>> long[]
		mapTypeSubscripts(final ImgPlus<T> hyperStack, final int[] splitIndices)
	{
		final long[] subscripts = new long[splitIndices.length];
		final Map<AxisType, Integer> typeCounts = new HashMap<>();
		for (int i = 0; i < splitIndices.length; i++) {
			final int splitIndex = splitIndices[i];
			final AxisType type = hyperStack.axis(splitIndex).type();
			typeCounts.compute(type, (key, count) -> count == null ? 1 : count + 1);
			subscripts[i] = typeCounts.get(type);
		}
		return subscripts;
	}

	/**
	 * Recursively calls {@link #applySplit(ImgPlus, List)} to split the
	 * hyperstack into subspaces
	 *
	 * @param hyperstack an n-dimensional image
	 * @param splitIndices the indices of the axes in the hyperstack used for
	 *          splitting
	 * @param splitIndex the i in splitIndices[i] currently used. Start from the
	 *          last index
	 * @param meta the metadata describing the position of the next subspace
	 * @param splitCoordinates the (dimension, position) pairs describing the
	 *          current split
	 * @param subscripts the subscripts of the axes see
	 * @param subspaces A builder for the stream of all the subspaces formed
	 */
	private static <T extends RealType<T> & NativeType<T>> void splitDims(
		final ImgPlus<T> hyperstack, final int[] splitIndices, final int splitIndex,
		final HyperAxisMeta[] meta,
		final List<ValuePair<IntType, LongType>> splitCoordinates,
		final long[] subscripts, final Builder<Subspace<T>> subspaces)
	{
		if (splitIndex < 0) {
			final RandomAccessibleInterval<T> subspace = applySplit(hyperstack,
				splitCoordinates);
			if (!isEmptySubspace(subspace)) {
				subspaces.add(new Subspace<>(subspace, meta));
			}
		}
		else {
			final int splitDimension = splitIndices[splitIndex];
			final AxisType type = hyperstack.axis(splitDimension).type();
			final long subscript = subscripts[splitIndex];
			final long size = hyperstack.dimension(splitDimension);
			final ValuePair<IntType, LongType> pair = new ValuePair<>(new IntType(
				splitDimension), new LongType());
			for (long position = 0; position < size; position++) {
				pair.b.set(position);
				splitCoordinates.add(pair);
				meta[splitIndex] = new HyperAxisMeta(type, position, subscript);
				splitDims(hyperstack, splitIndices, splitIndex - 1, meta,
					splitCoordinates, subscripts, subspaces);
				splitCoordinates.remove(pair);
			}
		}
	}

	private static <T extends RealType<T> & NativeType<T>> boolean
		isEmptySubspace(final RandomAccessibleInterval<T> hyperSlice)
	{
		return hyperSlice.numDimensions() == 0;
	}

	/**
	 * Splits a subspace along the given coordinates
	 * <p>
	 * For example, if you have a 5D {X, Y, Z, C, T} hyperstack, and give the
	 * coordinates {{3, 0}, {4, 1}} you'll get a 3D {X, Y, Z} subspace of the
	 * first channel, and second time frame
	 * </p>
	 * 
	 * @param hyperstack an n-dimensional image
	 * @param splitCoordinates (dimension, position) pairs describing the
	 *          hyperstack split
	 * @return The subspace interval
	 */
	private static <T extends RealType<T> & NativeType<T>>
		RandomAccessibleInterval<T> applySplit(final ImgPlus<T> hyperstack,
			final List<ValuePair<IntType, LongType>> splitCoordinates)
	{
		final List<ValuePair<IntType, LongType>> workingSplit = createWorkingCopy(
			splitCoordinates);
		RandomAccessibleInterval<T> slice = hyperstack;
		for (int i = 0; i < workingSplit.size(); i++) {
			final int dimension = workingSplit.get(i).a.get();
			final long position = workingSplit.get(i).b.get();
			slice = Views.hyperSlice(slice, dimension, position);
			decrementIndices(workingSplit, dimension);
		}
		return slice;
	}

	/**
	 * Clones and sorts the given {@link List}.
	 * <p>
	 * It ensures that original ones don't get altered while applying a split (see
	 * {@link #applySplit(ImgPlus, List)}). Pairs are sorted in the order of
	 * dimension ({@link ValuePair#a}).
	 * </p>
	 */
	private static List<ValuePair<IntType, LongType>> createWorkingCopy(
		final List<ValuePair<IntType, LongType>> splitCoordinates)
	{
		final List<ValuePair<IntType, LongType>> workingSplit = new ArrayList<>();
		for (ValuePair<IntType, LongType> pair : splitCoordinates) {
			ValuePair<IntType, LongType> copy = new ValuePair<>(pair.a.copy(),
				pair.b);
			workingSplit.add(copy);
		}
		workingSplit.sort(Comparator.comparingInt(pair -> pair.a.get()));
		return workingSplit;
	}

	/**
	 * A helper method for {@link #applySplit(ImgPlus, List)} that ensures that it
	 * doesn't throw a {@link IndexOutOfBoundsException}
	 * <p>
	 * After calling {@link Views#hyperSlice(RandomAccessibleInterval, int, long)}
	 * on an n-dimensional {@link ImgPlus} the resulting
	 * {@link RandomAccessibleInterval} will have n-1 dimensions. Thus the
	 * dimension indices in splitCoordinates need to be decremented if they come
	 * after the index used in the split.
	 * </p>
	 *
	 * @param splitCoordinates (dimension, position) pairs describing a hyperspace
	 *          split
	 * @param dimension The index of the dimension in the last split
	 */
	private static void decrementIndices(
		final List<ValuePair<IntType, LongType>> splitCoordinates,
		final int dimension)
	{
		for (ValuePair<IntType, LongType> pair : splitCoordinates) {
			IntType pairDimension = pair.getA();
			if (pairDimension.get() >= dimension) {
				pairDimension.dec();
			}
		}
	}

	/**
	 * Finds the axis used to split the hyperstack
	 * <p>
	 * For example, if you want to split a 5D {X, Y, C, Z, T} {@link ImgPlus} into
	 * 3D {X, Y, Z} subspaces, call with types = {Axes.X, Axes.Y, Axes.Z}
	 * </p>
	 *
	 * @param hyperstack An n-dimensional image
	 * @param types The axes that define the subspaces
	 * @return Indices of the axis used to split the hyperstack
	 */
	private static <T extends RealType<T> & NativeType<T>> int[]
		findSplitAxisIndices(ImgPlus<T> hyperstack, final List<AxisType> types)
	{
		final int n = hyperstack.numDimensions();
		return IntStream.range(0, n).filter(d -> !isAnyOfTypes(hyperstack.axis(d),
			types)).toArray();
	}

	private static boolean isAnyOfTypes(final TypedAxis axis,
		final List<AxisType> types)
	{
		return types.stream().anyMatch(t -> axis.type() == t);
	}
	// endregion

	// region -- Helper classes --

	/**
	 * A class which stores a subspace interval of an n-dimensional hyperspace,
	 * and metadata
	 * <p>
	 * The metadata describes the position of the subspace in the hyperspace.
	 * </p>
	 */
	public static final class Subspace<T extends RealType<T> & NativeType<T>> {

		public final RandomAccessibleInterval<T> interval;
		private final List<HyperAxisMeta> subspaceMeta;

		/**
		 * Creates a subspace record
		 *
		 * @param subspace An interval which defines a subspace
		 * @param subspaceMeta positions of the subspace in the hyperspace
		 */
		private Subspace(final RandomAccessibleInterval<T> subspace,
			final HyperAxisMeta[] subspaceMeta)
		{
			this.interval = subspace;
			if (subspaceMeta == null) {
				this.subspaceMeta = new ArrayList<>();
				return;
			}
			this.subspaceMeta = Arrays.stream(subspaceMeta).filter(Objects::nonNull)
				.collect(Collectors.toList());
		}

		/**
		 * Position of the subspace in each additional dimension of the hyperspace.
		 * <p>
		 * For example, one 3D {X, Y, Z} subspace of a 5D {X, Y, C, Z, T}
		 * hyperspace, would have position {0, 1} - 1st channel, 2nd frame.
		 * </p>
		 */
		public LongStream getPosition() {
			return subspaceMeta.stream().mapToLong(HyperAxisMeta::getPosition);
		}

		/**
		 * Types of the additional hyperspace dimensions
		 * <p>
		 * For example a 3D {X, Y, Z} subspace of a 5D {X, Y, C, Z, T} hyperspace,
		 * would have types {Axes.CHANNEL, Axes.TIME}.
		 * </p>
		 */
		public Stream<AxisType> getAxisTypes() {
			return subspaceMeta.stream().map(HyperAxisMeta::getType);
		}

		/**
		 * Subscripts of the additional hyperspace dimensions
		 * <p>
		 * Subscripts identify multiple axis of the same type
		 * </p>
		 * <p>
		 * For example a 3D {X, Y} subspace of a 6D {X, Y, C, Z, T, T} hyperspace,
		 * would have subscripts {1, 1, 1, 2}.
		 * </p>
		 */
		public LongStream getSubScripts() {
			return subspaceMeta.stream().mapToLong(HyperAxisMeta::getSubscript);
		}

		@Override
		public String toString() {
			return subspaceMeta.stream().map(HyperAxisMeta::toString).reduce((a,
				b) -> a + ", " + b).orElse("");
		}

		/**
		 * Describes the metadata of the subspace in relation to one of the axes in
		 * the hyperspace
		 */
		protected static final class HyperAxisMeta {

			/** {@link AxisType} of the dimension */
			private final AxisType type;
			/** The position of the subspace in the dimension */
			private final long position;
			/** An ID number to separate the dimension from others of the same type */
			private final long subscript;
			/**
			 * A number added to the position when it's printed in
			 * {@link #toString()}
			 */
			private final long stringOffset;

			private HyperAxisMeta(final AxisType type, final long position,
				final long subscript)
			{
				this.type = type;
				this.position = position;
				this.subscript = subscript;
				stringOffset = ResultUtils.toConventionalIndex(type, 0);
			}

			private AxisType getType() {
				return type;
			}

			private long getPosition() {
				return position;
			}

			private long getSubscript() {
				return subscript;
			}

			@Override
			public String toString() {
				final String typeIndex = subscript > 1 ? "(" + subscript + ")" : "";
				return type.toString() + typeIndex + ": " + (position + stringOffset);
			}
		}
	}
	// endregion
}
