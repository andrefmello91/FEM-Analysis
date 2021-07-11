using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Storage;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic quantity vector class.
	/// </summary>
	/// <typeparam name="TQuantity">The quantity that represents the value of components of the vector.</typeparam>
	/// <typeparam name="TUnit">The unit enumeration that represents the quantity of the components of this vector.</typeparam>
	public abstract class QuantityVector<TQuantity, TUnit> : DenseVector, IUnitConvertible<TUnit>, ICloneable<QuantityVector<TQuantity, TUnit>>, IEquatable<QuantityVector<TQuantity, TUnit>>, IEnumerable<TQuantity>
		where TQuantity : IQuantity<TUnit>
		where TUnit : Enum
	{

		#region Fields

		private TUnit _unit;

		#endregion

		#region Properties

		/// <summary>
		///     Get/set the quantity at this index.
		/// </summary>
		/// <param name="index">The index of the component.</param>
		public new TQuantity this[int index]
		{
			get => (TQuantity) Values[index].As(_unit);
			set => Values[index] = value.As(_unit);
		}

		#region Interface Implementations

		/// <inheritdoc />
		public TUnit Unit
		{
			get => _unit;
			set => ChangeUnit(value);
		}

		#endregion

		#endregion

		#region Constructors

		/// <inheritdoc />
		protected QuantityVector(IEnumerable<double> values, TUnit unit)
			: this(DenseVectorStorage<double>.OfEnumerable(values.ToArray()), unit)
		{
		}

		/// <inheritdoc />
		/// <param name="quantities">The component values of the vector.</param>
		protected QuantityVector(IEnumerable<TQuantity> quantities)
			: this(quantities.Select(v => v.As(quantities.First().Unit)), quantities.First().Unit)
		{
		}

		/// <inheritdoc />
		private QuantityVector(DenseVectorStorage<double> storage, TUnit unit)
			: base(storage) =>
			_unit = unit;

		#endregion

		#region Methods

		/// <inheritdoc cref="Vector{T}.AbsoluteMaximum" />
		public new TQuantity AbsoluteMaximum() => (TQuantity) base.AbsoluteMaximum().As(Unit);

		/// <inheritdoc cref="Vector{T}.AbsoluteMinimum" />
		public new TQuantity AbsoluteMinimum() => (TQuantity) base.AbsoluteMinimum().As(Unit);

		/// <inheritdoc cref="Vector{T}.Add(T)" />
		public QuantityVector<TQuantity, TUnit> Add(TQuantity scalar)
		{
			var result = Clone();

			if (scalar.Value.Equals(0))
				return result;

			result.Clear();
			DoAdd(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.Add(Vector{T})" />
		public QuantityVector<TQuantity, TUnit> Add(QuantityVector<TQuantity, TUnit> other)
		{
			var result = Count == other.Count
				? Clone()
				: throw new ArgumentException("All vectors must have the same dimensionality.", nameof(other));

			result.Clear();

			DoAdd(Unit.Equals(other.Unit)
					? other
					: other.Convert(Unit),
				result);

			return result;
		}

		/// <inheritdoc cref="IUnitConvertible{TUnit}.Convert" />
		public abstract QuantityVector<TQuantity, TUnit> Convert(TUnit unit);

		/// <inheritdoc cref="Vector{T}.Divide(T)" />
		public new QuantityVector<TQuantity, TUnit> Divide(double scalar)
		{
			var result = Clone();

			if (scalar.Equals(1))
				return result;

			result.Clear();

			DoDivide(scalar, result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.DivideByThis(T)" />
		public Vector<double> DivideByThis(TQuantity scalar)
		{
			var result = Build.Dense(Count);

			DoDivideByThis(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.DotProduct" />
		public double DotProduct(QuantityVector<TQuantity, TUnit> other) => Count == other.Count
			? DoDotProduct(other.Unit.Equals(Unit)
				? other
				: other.Convert(Unit))
			: throw new ArgumentException("All vectors must have the same dimensionality.", nameof(other));

		/// <inheritdoc cref="Vector{T}.Maximum" />
		public TQuantity Maximum() => (TQuantity) base.Maximum().As(Unit);

		/// <inheritdoc cref="Vector{T}.Minimum" />
		public TQuantity Minimum() => (TQuantity) base.Minimum().As(Unit);

		/// <inheritdoc cref="Vector{T}.Multiply(T)" />
		public new QuantityVector<TQuantity, TUnit> Multiply(double scalar)
		{
			var result = Clone();

			if (scalar.Equals(1))
				return result;

			result.Clear();

			if (scalar.Equals(0))
				return result;

			DoMultiply(scalar, result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.Negate()" />
		public new QuantityVector<TQuantity, TUnit> Negate()
		{
			var result = Clone();

			result.Clear();

			DoNegate(result);

			return result;
		}

		/// <summary>
		///     Get the vector simplified by constrained DoFs.
		/// </summary>
		/// <remarks>
		///     This uses the default tolerance.
		/// </remarks>
		/// <returns>
		///     The simplified <see cref="Vector{T}" />.
		/// </returns>
		public abstract QuantityVector<TQuantity, TUnit> Simplified(IEnumerable<int>? indexes = null);

		/// <summary>
		///     Get the vector simplified by constrained DoFs.
		/// </summary>
		/// <param name="threshold">
		///     A value for setting all values whose absolute value is smaller than to zero. If null, this is
		///     not applied.
		/// </param>
		/// <param name="indexes">The collection of indexes to simplify.</param>
		/// <returns>
		///     The simplified <see cref="Vector{T}" />.
		/// </returns>
		public abstract QuantityVector<TQuantity, TUnit> Simplified(double? threshold, IEnumerable<int>? indexes = null);

		/// <inheritdoc cref="Simplified(double?, IEnumerable{int}?)" />
		public QuantityVector<TQuantity, TUnit> Simplified(TQuantity? threshold, IEnumerable<int>? indexes = null) => Simplified(threshold?.As(Unit), indexes);

		/// <inheritdoc cref="Vector{T}.Subtract(T)" />
		public QuantityVector<TQuantity, TUnit> Subtract(TQuantity scalar)
		{
			var result = Clone();

			if (scalar.Value.Equals(0))
				return result;

			result.Clear();
			DoSubtract(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.Subtract(Vector{T})" />
		public QuantityVector<TQuantity, TUnit> Subtract(QuantityVector<TQuantity, TUnit> other)
		{
			var result = Count == other.Count
				? Clone()
				: throw new ArgumentException("All vectors must have the same dimensionality.", nameof(other));

			result.Clear();

			DoSubtract(Unit.Equals(other.Unit)
					? other
					: other.Convert(Unit),
				result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.SubtractFrom(T)" />
		public QuantityVector<TQuantity, TUnit> SubtractFrom(TQuantity scalar)
		{
			var result = Clone();

			if (scalar.Value.Equals(0))
				return result;

			result.Clear();
			DoSubtractFrom(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Vector{T}.Sum" />
		public new TQuantity Sum() => (TQuantity) base.Sum().As(Unit);

		#region Interface Implementations

		/// <inheritdoc cref="ICloneable{T}.Clone" />
		public abstract QuantityVector<TQuantity, TUnit> Clone();

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

		/// <inheritdoc />
		public IEnumerator<TQuantity> GetEnumerator() => Values
			.GetQuantities<TQuantity, TUnit>(Unit)
			.GetEnumerator();

		/// <inheritdoc />
		public bool Equals(QuantityVector<TQuantity, TUnit>? other)
		{
			if (other is null)
				return false;

			return
				base.Equals(Unit.Equals(other.Unit)
					? other
					: other.Convert(Unit));
		}

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;

			for (var i = 0; i < Count; i++)
				Values[i] = this[i].As(unit);

			_unit = unit;
		}

		IUnitConvertible<TUnit> IUnitConvertible<TUnit>.Convert(TUnit unit) => Convert(unit);

		#endregion

		#region Object override

		/// <inheritdoc cref="object.Equals(object)" />
		public new bool Equals(object? obj) =>
			obj is QuantityVector<TQuantity, TUnit> other && Equals(other);

		/// <inheritdoc cref="Add(TQuantity)" />
		public static QuantityVector<TQuantity, TUnit> operator +(QuantityVector<TQuantity, TUnit> left, TQuantity right) => left.Add(right);

		/// <inheritdoc cref="Add(QuantityVector{TQuantity, TUnit})" />
		public static QuantityVector<TQuantity, TUnit> operator +(QuantityVector<TQuantity, TUnit> left, QuantityVector<TQuantity, TUnit> right) => left.Add(right);

		/// <inheritdoc cref="Divide(double)" />
		public static QuantityVector<TQuantity, TUnit> operator /(QuantityVector<TQuantity, TUnit> left, double divisor) => left.Divide(divisor);

		/// <inheritdoc cref="DivideByThis(TQuantity)" />
		public static Vector<double> operator /(TQuantity left, QuantityVector<TQuantity, TUnit> divisor) => divisor.DivideByThis(left);

		/// <inheritdoc cref="StiffnessMatrix.op_Equality" />
		public static bool operator ==(QuantityVector<TQuantity, TUnit>? left, QuantityVector<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <inheritdoc cref="StiffnessMatrix.op_Inequality" />
		public static bool operator !=(QuantityVector<TQuantity, TUnit>? left, QuantityVector<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		/// <inheritdoc cref="Multiply(double)" />
		public static QuantityVector<TQuantity, TUnit> operator *(QuantityVector<TQuantity, TUnit> left, double multiplier) => left.Multiply(multiplier);

		/// <inheritdoc cref="Multiply(double)" />
		public static QuantityVector<TQuantity, TUnit> operator *(double multiplier, QuantityVector<TQuantity, TUnit> right) => right.Multiply(multiplier);

		/// <returns>
		///     The dot product between the vectors.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static double operator *(QuantityVector<TQuantity, TUnit> left, QuantityVector<TQuantity, TUnit> right) => left.DotProduct(right);

		/// <inheritdoc cref="Subtract(TQuantity)" />
		public static QuantityVector<TQuantity, TUnit> operator -(QuantityVector<TQuantity, TUnit> left, TQuantity right) => left.Subtract(right);

		/// <inheritdoc cref="Subtract(QuantityVector{TQuantity, TUnit})" />
		public static QuantityVector<TQuantity, TUnit> operator -(QuantityVector<TQuantity, TUnit> left, QuantityVector<TQuantity, TUnit> right) => left.Subtract(right);

		/// <inheritdoc cref="SubtractFrom(TQuantity)" />
		public static QuantityVector<TQuantity, TUnit> operator -(TQuantity left, QuantityVector<TQuantity, TUnit> right) => right.SubtractFrom(left);

		/// <inheritdoc cref="Negate()" />
		public static QuantityVector<TQuantity, TUnit> operator -(QuantityVector<TQuantity, TUnit> vector) => vector.Negate();

		/// <inheritdoc cref="object.ToString" />
		public new string ToString() =>
			$"Unit: {Unit} \n" +
			$"Value: {base.ToString()}";

		#endregion

		#endregion

	}
}