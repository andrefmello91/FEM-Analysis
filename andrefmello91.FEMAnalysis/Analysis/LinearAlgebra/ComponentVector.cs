using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic component vector class.
	/// </summary>
	/// <typeparam name="TQuantity">The quantity that represents the value of components of the vector.</typeparam>
	/// <typeparam name="TUnit">The unit enumeration that represents the quantity of the components of this vector.</typeparam>
	public abstract class ComponentVector<TQuantity, TUnit> : IUnitConvertible<TUnit>, ICloneable<ComponentVector<TQuantity, TUnit>>, IEquatable<ComponentVector<TQuantity, TUnit>>, IEnumerable<TQuantity>
		where TQuantity : IQuantity<TUnit>
		where TUnit : Enum
	{

		#region Fields

		/// <summary>
		///     Get/set the matrix value of this object, with components in <see cref="Unit" />.
		/// </summary>
		protected readonly List<double> Values;

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

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}(IEnumerable{TQuantity})" />
		/// <param name="unit">The unit of <paramref name="values" />'s components.</param>
		protected ComponentVector(IEnumerable<double> values, TUnit unit)
		{
			Values = new List<double>(values);
			_unit  = unit;
		}

		/// <summary>
		///     Create a component vector.
		/// </summary>
		/// <param name="values">The enumerable of vector's components.</param>
		protected ComponentVector(IEnumerable<TQuantity> values)
		{
			_unit = values.First().Unit;

			Values = values
				.Select(v => v.As(_unit))
				.ToList();
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="Vector{T}.AbsoluteMaximum" />
		public TQuantity AbsoluteMaximum() => (TQuantity) Values.MaximumAbsolute().As(Unit);

		/// <inheritdoc cref="Vector{T}.AbsoluteMinimum" />
		public TQuantity AbsoluteMinimum() => (TQuantity) Values.MinimumAbsolute().As(Unit);

		/// <inheritdoc cref="Vector{T}.Clear" />
		public void Clear()
		{
			for (var i = 0; i < Count; i++)
				Values[i] = 0;
		}

		/// <inheritdoc cref="IUnitConvertible{TUnit}.Convert" />
		public abstract ComponentVector<TQuantity, TUnit> Convert(TUnit unit);

		/// <inheritdoc cref="Vector{T}.Maximum" />
		public TQuantity Maximum() => (TQuantity) Values.Maximum().As(Unit);

		/// <inheritdoc cref="Vector{T}.Minimum" />
		public TQuantity Minimum() => (TQuantity) Values.Minimum().As(Unit);

		/// <inheritdoc cref="Vector{T}.Norm" />
		public double Norm(double p) => Values.ToVector().Norm(p);

		/// <summary>
		///     Get the vector simplified by constrained DoFs.
		/// </summary>
		/// <remarks>
		///		This uses the default tolerance.
		/// </remarks>
		/// <returns>
		///     The simplified <see cref="Vector{T}" />.
		/// </returns>
		public abstract Vector<double> Simplified();

		/// <summary>
		///     Get the vector simplified by constrained DoFs.
		/// </summary>
		/// <param name="threshold">
		///     A value for setting all values whose absolute value is smaller than to zero. If null, this is
		///     not applied.
		/// </param>
		/// <returns>
		///     The simplified <see cref="Vector{T}" />.
		/// </returns>
		public Vector<double> Simplified(double? threshold)
		{
			var simplified = Values.ToVector();

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
		public Matrix<double> ToColumnMatrix() => Values.ToVector().ToColumnMatrix();

		/// <inheritdoc cref="Vector{T}.ToRowMatrix" />
		public Matrix<double> ToRowMatrix() => Values.ToVector().ToRowMatrix();

		#region Interface Implementations

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;

			for (var i = 0; i < Count; i++)
				Values[i] = this[i].As(unit);

			_unit = unit;
		}

		/// <inheritdoc cref="ICloneable{T}.Clone" />
		public abstract ComponentVector<TQuantity, TUnit> Clone();

		/// <inheritdoc />
		public bool Equals(ComponentVector<TQuantity, TUnit>? other) =>
			other is not null && Values.ToVector() == other.Convert(Unit).Values.ToVector();

		/// <inheritdoc />
		public IEnumerator<TQuantity> GetEnumerator() => Values
			.GetQuantities<TQuantity, TUnit>(Unit)
			.GetEnumerator();

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
			$"Value: {Values.ToVector()}";

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
		///     The dot product between the vectors.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static double operator *(ComponentVector<TQuantity, TUnit> left, ComponentVector<TQuantity, TUnit> right) => left.ToVector(left.Unit) * right.ToVector(left.Unit);

		#endregion

	}
}