using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic component vector class.
	/// </summary>
	/// <typeparam name="TQuantity">The quantity that represents the value of components of the vector.</typeparam>
	/// <typeparam name="TUnit">The unit enumeration that represents the quantity of the components of this vector.</typeparam>
	public class ComponentVector<TQuantity, TUnit> : IUnitConvertible<TUnit>, ICloneable<ComponentVector<TQuantity, TUnit>>, IEquatable<ComponentVector<TQuantity, TUnit>>, IEnumerable<TQuantity>
		where TQuantity : IQuantity<TUnit>
		where TUnit : Enum
	{

		#region Fields

		/// <summary>
		///     Get/set the matrix value of this object, with components in <see cref="Unit" />.
		/// </summary>
		protected readonly List<TQuantity> Values;

		private TUnit _unit;

		#endregion

		#region Properties

		/// <summary>
		///     The index of constrained DoFs.
		/// </summary>
		public List<int>? ConstraintIndex { get; set; }

		/// <inheritdoc cref="Vector{T}.Count" />
		public int Count => Values.Count;

		/// <summary>
		///     Get/set the quantity at this index.
		/// </summary>
		/// <param name="index">The index of the component.</param>
		public TQuantity this[int index]
		{
			get => (TQuantity) Values[index].ToUnit(_unit);
			set => Values[index] = (TQuantity) value.ToUnit(_unit);
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

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}(IEnumerable{TQuantity})" />
		/// <param name="unit">The unit of <paramref name="values" />'s components.</param>
		public ComponentVector(IEnumerable<double> values, TUnit unit)
		{
			Values = values
				.GetQuantities<TQuantity, TUnit>(unit)
				.ToList();

			_unit = unit;
		}

		/// <summary>
		///     Create a component vector.
		/// </summary>
		/// <param name="values">The enumerable of vector's components.</param>
		public ComponentVector(IEnumerable<TQuantity> values)
		{
			_unit = values.First().Unit;

			Values = values
				.Select(v => v.ToUnit(_unit))
				.Cast<TQuantity>()
				.ToList();
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="Vector{T}.AbsoluteMaximum" />
		public TQuantity AbsoluteMaximum() => Values
			.Select(v => v.Abs())
			.Max(Unit);

		/// <inheritdoc cref="Vector{T}.AbsoluteMinimum" />
		public TQuantity AbsoluteMinimum() => Values
			.Select(v => v.Abs())
			.Min(Unit);

		/// <inheritdoc cref="Vector{T}.Clear" />
		public void Clear() => Values.Clear();

		/// <inheritdoc cref="IUnitConvertible{TUnit}.Convert" />
		public ComponentVector<TQuantity, TUnit> Convert(TUnit unit) => new(Values.Select(v => v.As(unit)), unit)
		{
			ConstraintIndex = ConstraintIndex
		};

		/// <inheritdoc cref="Vector{T}.Maximum" />
		public TQuantity Maximum() => Values.Max(Unit);

		/// <inheritdoc cref="Vector{T}.Minimum" />
		public TQuantity Minimum() => Values.Min(Unit);

		/// <inheritdoc cref="Vector{T}.Norm" />
		public double Norm(double p) => Values.ToVector(Unit).Norm(p);


		/// <summary>
		///     Get the vector simplified by constraint indexes.
		/// </summary>
		/// <param name="threshold">
		///     A value for setting all values whose absolute value is smaller than to zero. If null, this is
		///     not applied.
		/// </param>
		/// <returns>
		///     The simplified <see cref="Vector{T}" />.
		/// </returns>
		public Vector<double> Simplified(double? threshold = null)
		{
			var simplified = Values.ToVector(Unit);

			if (ConstraintIndex is not null)
				foreach (var index in ConstraintIndex)
					simplified[index] = 0;

			if (threshold.HasValue)
				simplified.CoerceZero(threshold.Value);

			return simplified;
		}

		/// <inheritdoc cref="Simplified(double?)" />
		public Vector<double> Simplified(TQuantity? threshold) => Simplified(threshold?.As(Unit));

		/// <inheritdoc cref="Vector{T}.ToColumnMatrix" />
		public Matrix<double> ToColumnMatrix() => Values.ToVector(Unit).ToColumnMatrix();

		/// <inheritdoc cref="Vector{T}.ToRowMatrix" />
		public Matrix<double> ToRowMatrix() => Values.ToVector(Unit).ToRowMatrix();

		#region Interface Implementations

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;

			for (var i = 0; i < Count; i++)
				Values[i] = (TQuantity) Values[i].ToUnit(unit);

			_unit = unit;
		}

		/// <inheritdoc cref="ICloneable{T}.Clone" />
		public ComponentVector<TQuantity, TUnit> Clone() => new(Values)
		{
			ConstraintIndex = ConstraintIndex
		};

		/// <inheritdoc />
		public bool Equals(ComponentVector<TQuantity, TUnit>? other) =>
			other is not null && _unit.Equals(other._unit) && Values.ToVector(Unit) == other.Values.ToVector(Unit);

		/// <inheritdoc />
		public IEnumerator<TQuantity> GetEnumerator() => Values.GetEnumerator();

		IUnitConvertible<TUnit> IUnitConvertible<TUnit>.Convert(TUnit unit) => Convert(unit);

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

		#endregion

		#region Object override

		/// <inheritdoc />
		public override bool Equals(object? obj) =>
			obj is ComponentVector<TQuantity, TUnit> other && Equals(other);

		/// <inheritdoc />
		public override int GetHashCode() => _unit.GetHashCode() * Values.GetHashCode();

		/// <inheritdoc />
		public override string ToString() =>
			$"Unit: {Unit} \n" +
			$"Value: {Values}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the corresponding <see cref="Vector{T}" />.
		/// </summary>
		public static implicit operator Vector<double>(ComponentVector<TQuantity, TUnit> vector) => vector.ToVector(vector.Unit);

		/// <inheritdoc cref="StiffnessMatrix.op_Equality" />
		public static bool operator ==(ComponentVector<TQuantity, TUnit>? left, ComponentVector<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <inheritdoc cref="StiffnessMatrix.op_Inequality" />
		public static bool operator !=(ComponentVector<TQuantity, TUnit>? left, ComponentVector<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		/// <returns>
		///     A new vector with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ComponentVector<TQuantity, TUnit> operator +(ComponentVector<TQuantity, TUnit> left, ComponentVector<TQuantity, TUnit> right) =>
			new(left.ToVector(left.Unit) + right.ToVector(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};


		/// <returns>
		///     A new vector with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ComponentVector<TQuantity, TUnit> operator -(ComponentVector<TQuantity, TUnit> left, ComponentVector<TQuantity, TUnit> right) =>
			new(left.ToVector(left.Unit) - right.ToVector(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};


		/// <returns>
		///     A vector with components multiplied by a value
		/// </returns>
		public static ComponentVector<TQuantity, TUnit> operator *(double multiplier, ComponentVector<TQuantity, TUnit> vector) =>
			new(vector.Values.Select(v => v.Value * multiplier), vector.Unit)
			{
				ConstraintIndex = vector.ConstraintIndex
			};


		/// <inheritdoc cref="op_Multiply(double, ComponentVector{TQuantity,TUnit}) " />
		public static ComponentVector<TQuantity, TUnit> operator *(ComponentVector<TQuantity, TUnit> vector, double multiplier) => multiplier * vector;

		/// <returns>
		///     The dot product between the vectors.
		/// </returns>
		/// <inheritdoc cref="op_Subtraction" />
		public static double operator *(ComponentVector<TQuantity, TUnit> left, ComponentVector<TQuantity, TUnit> right) => left.ToVector(left.Unit) * right.ToVector(left.Unit);

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static ComponentVector<TQuantity, TUnit> operator -(ComponentVector<TQuantity, TUnit> vector) => new(vector.Select(v => -v.Value), vector.Unit)
		{
			ConstraintIndex = vector.ConstraintIndex
		};


		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static ComponentVector<TQuantity, TUnit> operator /(ComponentVector<TQuantity, TUnit> vector, double divisor) => new(vector.Select(v => v.Value / divisor), vector.Unit)
		{
			ConstraintIndex = vector.ConstraintIndex
		};

		#endregion

	}
}