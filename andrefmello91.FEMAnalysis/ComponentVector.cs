﻿using System;
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

		private TUnit _unit;
		protected Vector<double> Value;

		#endregion

		#region Properties

		/// <summary>
		///     The index of constrained DoFs.
		/// </summary>
		public int[]? ConstraintIndex { get; set; }

		/// <inheritdoc cref="Vector{T}.Count" />
		public int Count => Value.Count;

		/// <summary>
		///     Get/set the quantity at this index.
		/// </summary>
		/// <param name="index">The index of the component.</param>
		public TQuantity this[int index]
		{
			get => (TQuantity) Value[index].As(_unit);
			set => Value[index] = value.As(_unit);
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

		/// <summary>
		///     Create a component vector.
		/// </summary>
		/// <param name="value">The double enumerable.</param>
		/// <param name="unit">The unit of <paramref name="value" />'s components.</param>
		public ComponentVector(IEnumerable<double> value, TUnit unit)
		{
			Value = value is Vector<double> vector
				? vector
				: value.ToVector();

			_unit = unit;
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="IUnitConvertible{TUnit}.Convert" />
		public ComponentVector<TQuantity, TUnit> Convert(TUnit unit) => new(Value * Quantity.From(1, _unit).As(unit), unit);

		#region Interface Implementations

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;

			Value *= Quantity.From(1, _unit).As(unit);
			_unit =  unit;
		}

		/// <inheritdoc cref="ICloneable{T}.Clone" />
		public ComponentVector<TQuantity, TUnit> Clone() => new(Value.Clone(), _unit);

		/// <inheritdoc />
		public bool Equals(ComponentVector<TQuantity, TUnit>? other) =>
			other is not null && _unit.Equals(other._unit) && Value == other.Value;

		/// <inheritdoc />
		public IEnumerator<TQuantity> GetEnumerator() => Value
			.Select(v => (TQuantity) v.As(_unit))
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
		public override int GetHashCode() => _unit.GetHashCode() * Value.GetHashCode();

		/// <inheritdoc />
		public override string ToString() =>
			$"Unit: {Unit} \n" +
			$"Value: {Value}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the corresponding <see cref="Vector{T}" />.
		/// </summary>
		public static implicit operator Vector<double>(ComponentVector<TQuantity, TUnit> vector) => vector.Value;

		/// <inheritdoc cref="StiffnessMatrix.op_Equality" />
		public static bool operator ==(ComponentVector<TQuantity, TUnit>? left, ComponentVector<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <inheritdoc cref="StiffnessMatrix.op_Inequality" />
		public static bool operator !=(ComponentVector<TQuantity, TUnit>? left, ComponentVector<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		/// <returns>
		///     A new vector with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ComponentVector<TQuantity, TUnit> operator +(ComponentVector<TQuantity, TUnit> left, ComponentVector<TQuantity, TUnit> right) => new(left.Value + right.Convert(left.Unit).Value, left.Unit);

		/// <returns>
		///     A new vector with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ComponentVector<TQuantity, TUnit> operator -(ComponentVector<TQuantity, TUnit> left, ComponentVector<TQuantity, TUnit> right) => new(left.Value - right.Convert(left.Unit).Value, left.Unit);

		/// <returns>
		///     A vector with components multiplied by a value
		/// </returns>
		public static ComponentVector<TQuantity, TUnit> operator *(double value, ComponentVector<TQuantity, TUnit> right) => new(value * right.Value, right.Unit);

		/// <inheritdoc cref="op_Multiply(double, ComponentVector{TQuantity,TUnit}) " />
		public static ComponentVector<TQuantity, TUnit> operator *(ComponentVector<TQuantity, TUnit> left, double value) => value * left;

		#endregion

	}
}