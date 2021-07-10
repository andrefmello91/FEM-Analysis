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

		/// <inheritdoc/>
		protected QuantityVector(IEnumerable<double> values, TUnit unit)
			: this(DenseVectorStorage<double>.OfEnumerable(values), unit)
		{
		}
		
		/// <inheritdoc/>
		/// <param name="quantities">The component values of the vector.</param>
		protected QuantityVector(IEnumerable<TQuantity> quantities)
			: this(quantities.Select(v => v.As(quantities.First().Unit)), quantities.First().Unit)
		{
		}
		
		/// <inheritdoc />
		private QuantityVector(DenseVectorStorage<double> storage, TUnit unit)
			: base(storage)
		{
			_unit = unit;
		}
		
		#endregion

		#region Methods

		/// <inheritdoc cref="Vector{T}.AbsoluteMaximum" />
		public new TQuantity AbsoluteMaximum() => (TQuantity) base.AbsoluteMaximum().As(Unit);

		/// <inheritdoc cref="Vector{T}.AbsoluteMinimum" />
		public new TQuantity AbsoluteMinimum() => (TQuantity) base.AbsoluteMinimum().As(Unit);

		/// <inheritdoc cref="IUnitConvertible{TUnit}.Convert" />
		public abstract QuantityVector<TQuantity, TUnit> Convert(TUnit unit);

		/// <inheritdoc cref="Vector{T}.Maximum" />
		public TQuantity Maximum() => (TQuantity) base.Maximum().As(Unit);

		/// <inheritdoc cref="Vector{T}.Minimum" />
		public TQuantity Minimum() => (TQuantity) base.Minimum().As(Unit);

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
			
			other = Unit.Equals(other.Unit)
				? other
				: other.Convert(Unit);

			return
				base.Equals(other);
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

		/// <inheritdoc cref="object.Equals(object)"/>
		public new bool Equals(object? obj) =>
			obj is QuantityVector<TQuantity, TUnit> other && Equals(other);

		/// <inheritdoc cref="StiffnessMatrix.op_Equality" />
		public static bool operator ==(QuantityVector<TQuantity, TUnit>? left, QuantityVector<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <inheritdoc cref="StiffnessMatrix.op_Inequality" />
		public static bool operator !=(QuantityVector<TQuantity, TUnit>? left, QuantityVector<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		/// <returns>
		///     The dot product between the vectors.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static double operator *(QuantityVector<TQuantity, TUnit> left, QuantityVector<TQuantity, TUnit> right)
		{
			var other = right.Unit.Equals(left.Unit)
				? right
				: right.Convert(left.Unit);

			return
				(Vector<double>) left * other;
		}

		/// <inheritdoc cref="object.ToString" />
		public new string ToString() =>
			$"Unit: {Unit} \n" +
			$"Value: {base.ToString()}";

		#endregion

		#endregion

	}
}