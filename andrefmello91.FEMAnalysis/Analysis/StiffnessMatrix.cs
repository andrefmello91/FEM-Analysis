using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Generic stiffness matrix class.
	/// </summary>
	/// <typeparam name="TQuantity">The quantity that represents the value of components of the matrix.</typeparam>
	/// <typeparam name="TUnit">The unit enumeration that represents the quantity of the components of the matrix.</typeparam>
	public class StiffnessMatrix<TQuantity, TUnit> : IUnitConvertible<TUnit>, ICloneable<StiffnessMatrix<TQuantity, TUnit>>, IEquatable<StiffnessMatrix<TQuantity, TUnit>>
		where TQuantity : IQuantity<TUnit>
		where TUnit : Enum
	{
		
		#region Fields

		private TUnit _unit;

		/// <summary>
		///     The corresponding matrix, with components in <see cref="Unit" />.
		/// </summary>
		protected Matrix<double> Value;

		#endregion

		#region Properties

		#region Interface Implementations

		/// <inheritdoc />
		public TUnit Unit
		{
			get => _unit;
			set => ChangeUnit(value);
		}

		/// <summary>
		///		The index of constrained DoFs.
		/// </summary>
		public int[]? ConstraintIndex { get; set; }

		/// <inheritdoc cref="Matrix{T}.RowCount"/>
		public int Rows => Value.RowCount;
		
		/// <inheritdoc cref="Matrix{T}.ColumnCount"/>
		public int Columns => Value.ColumnCount;

		/// <summary>
		///		Get/set the value at these indexes.
		/// </summary>
		/// <param name="rowIndex">The row of the required element.</param>
		/// <param name="columnIndex">The column of the required element.</param>
		public TQuantity this[int rowIndex, int columnIndex]
		{
			get => (TQuantity) Value[rowIndex, columnIndex].As(_unit);
			set => Value[rowIndex, columnIndex] = value.As(_unit);
		}

		#endregion

		#endregion

		#region Constructors

		/// <inheritdoc cref="StiffnessMatrix{T,T}(Matrix{double}, TUnit)" />
		public StiffnessMatrix(double[,] value, TUnit unit)
			: this(Matrix<double>.Build.DenseOfArray(value), unit)
		{
		}

		/// <summary>
		///     Create a stiffness matrix.
		/// </summary>
		/// <param name="value">The <see cref="Matrix{T}" /> or <see cref="double" /> array value.</param>
		/// <param name="unit">The unit of <paramref name="value" />'s components</param>
		public StiffnessMatrix(Matrix<double> value, TUnit unit)
		{
			Value = value;
			_unit  = unit;
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public StiffnessMatrix<TQuantity, TUnit> Convert(TUnit unit) => new(Quantity.From(1, Unit).As(unit) * Value, unit);

		#region Interface Implementations

		/// <inheritdoc cref="Matrix{T}.Row(int)"/>
		public Vector<double> Row(int index) => Value.Row(index);
		
		/// <inheritdoc cref="Matrix{T}.Column(int)"/>
		public Vector<double> Column(int index) => Value.Column(index);

		/// <inheritdoc cref="Matrix{T}.ClearRows(int[])"/>
		public void ClearRows(params int[] indexes) => Value.ClearRows(indexes);
		
		/// <inheritdoc cref="Matrix{T}.ClearColumns(int[])"/>
		public void ClearColumns(params int[] indexes) => Value.ClearColumns(indexes);

		/// <inheritdoc cref="Matrix{T}.Clear"/>
		public void Clear() => Value.Clear();
		
		/// <inheritdoc cref="Matrix{T}.Determinant"/>
		public double Determinant() => Value.Determinant();
		
		/// <inheritdoc cref="Matrix{T}.Transpose()"/>
		public StiffnessMatrix<TQuantity, TUnit> Transpose() => new (Value.Transpose(), _unit);

		/// <summary>
		///		Get the simplified stiffness matrix by the constrained DoFs.
		/// </summary>
		public Matrix<double> Simplified() => ConstraintIndex is not null
			? StiffnessMatrix.SimplifiedStiffness(Value, ConstraintIndex)
			: Value;
		
		/// <summary>
		///		Transform this stiffness to another coordinate system.
		/// </summary>
		/// <param name="transformationMatrix">The transformation matrix.</param>
		/// <returns>
		///		A new <see cref="StiffnessMatrix"/> with transformed components.
		/// </returns>
		/// <exception cref="ArgumentException">If the dimensions of <paramref name="transformationMatrix"/> don't conform with this.</exception>
		public StiffnessMatrix<TQuantity, TUnit> Transform(Matrix<double> transformationMatrix)
		{
			var value = transformationMatrix.Transpose() * Value * transformationMatrix;
			
			return
				new StiffnessMatrix<TQuantity, TUnit>(value, _unit);
		}

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;
			
			// Multiply matrix
			Value *= Quantity.From(1, Unit).As(unit);
			
			// Set
			_unit = unit;
		}

		/// <inheritdoc />
		public StiffnessMatrix<TQuantity, TUnit> Clone() => new(Value.Clone(), _unit);

		/// <inheritdoc />
		public bool Equals(StiffnessMatrix<TQuantity, TUnit>? other) =>
			other is not null && _unit.Equals(other._unit) && Value.Equals(other.Value);

		/// <inheritdoc />
		IUnitConvertible<TUnit> IUnitConvertible<TUnit>.Convert(TUnit unit) => Convert(unit);

		#endregion

		#region Object override

		/// <inheritdoc />
		public override bool Equals(object? obj) =>
			obj is StiffnessMatrix<TQuantity, TUnit> other && Equals(other);

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
		///     Get the corresponding <see cref="Matrix{T}" /> value of a <see cref="StiffnessMatrix" />.
		/// </summary>
		public static implicit operator Matrix<double>(StiffnessMatrix<TQuantity, TUnit> stiffnessMatrix) => stiffnessMatrix.Value;

		/// <returns>
		///     True if objects are equal.
		/// </returns>
		public static bool operator ==(StiffnessMatrix<TQuantity, TUnit>? left, StiffnessMatrix<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <returns>
		///     True if objects are not equal.
		/// </returns>
		public static bool operator !=(StiffnessMatrix<TQuantity, TUnit>? left, StiffnessMatrix<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		/// <returns>
		///     A new stiffness matrix with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentOutOfRangeException">If left and right don't have the same dimensions.</exception>
		public static StiffnessMatrix<TQuantity, TUnit> operator +(StiffnessMatrix<TQuantity, TUnit> left, StiffnessMatrix<TQuantity, TUnit> right) => new(left.Value + right.Convert(left.Unit).Value, left.Unit);

		/// <returns>
		///     A new stiffness matrix with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentOutOfRangeException">If left and right don't have the same dimensions.</exception>
		public static StiffnessMatrix<TQuantity, TUnit> operator -(StiffnessMatrix<TQuantity, TUnit> left, StiffnessMatrix<TQuantity, TUnit> right) => new(left.Value - right.Convert(left.Unit).Value, left.Unit);

		/// <returns>
		///     A new stiffness matrix with components multiplied by a value
		/// </returns>
		public static StiffnessMatrix<TQuantity, TUnit> operator *(double value, StiffnessMatrix<TQuantity, TUnit> right) => new(value * right.Value, right.Unit);

		/// <inheritdoc cref="op_Multiply(double, StiffnessMatrix{TQuantity,TUnit}) " />
		public static StiffnessMatrix<TQuantity, TUnit> operator *(StiffnessMatrix<TQuantity, TUnit> left, double value) => value * left;

		#endregion
	}
	
	/// <summary>
	///     Default stiffness matrix class.
	/// </summary>
	/// <remarks>
	///		Unit is <see cref="ForcePerLengthUnit.NewtonPerMillimeter"/>.
	///		<para>
	///		Quantity is <see cref="ForcePerLength"/>.
	///		</para>
	/// </remarks>
	public class StiffnessMatrix : StiffnessMatrix<ForcePerLength, ForcePerLengthUnit>
	{
		#region Constructors

		/// <inheritdoc cref="StiffnessMatrix(Matrix{double}, ForcePerLengthUnit)" />
		public StiffnessMatrix(double[,] value, ForcePerLengthUnit unit = ForcePerLengthUnit.NewtonPerMillimeter)
			: this(Matrix<double>.Build.DenseOfArray(value), unit)
		{
		}

		/// <summary>
		///     Create a stiffness matrix.
		/// </summary>
		/// <param name="value">The <see cref="Matrix{T}" /> or <see cref="double" /> array value.</param>
		/// <param name="unit">The unit of <paramref name="value" />'s components</param>
		public StiffnessMatrix(Matrix<double> value, ForcePerLengthUnit unit = ForcePerLengthUnit.NewtonPerMillimeter)
			: base(value, unit)
		{
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public new StiffnessMatrix Convert(ForcePerLengthUnit unit) => (StiffnessMatrix) base.Convert(unit);

		/// <inheritdoc cref="Matrix{T}.Transpose()"/>
		public new StiffnessMatrix Transpose() => (StiffnessMatrix) base.Transpose();

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.Transform"/>
		public new StiffnessMatrix Transform(Matrix<double> transformationMatrix) => (StiffnessMatrix) base.Transform(transformationMatrix);

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.Clone"/>
		public new StiffnessMatrix Clone() => (StiffnessMatrix) base.Clone();

		#endregion

		#region Object override

		#endregion

		#region Operators

		/// <summary>
		///     Create a <see cref="StiffnessMatrix" /> from a <see cref="Matrix{T}" />, with components in
		///     <see cref="ForcePerLengthUnit.NewtonPerMillimeter" />.
		/// </summary>
		public static implicit operator StiffnessMatrix(Matrix<double> value) => new(value);

		#endregion

		/// <summary>
		///     Get the global stiffness simplified.
		/// </summary>
		/// <param name="stiffness">The global stiffness <see cref="Matrix{T}" /> to simplify.</param>
		/// <param name="indexes">The DoF indexes to simplify matrix.</param>
		/// <param name="simplifyZeroRows">Simplify matrix at rows containing only zero elements?</param>
		public static Matrix<double> SimplifiedStiffness(Matrix<double> stiffness, IEnumerable<int> indexes, bool simplifyZeroRows = true)
		{
			var simplifiedStiffness = stiffness.Clone();

			var index = indexes.ToArray();

			// Clear the rows and columns in the stiffness matrix
			simplifiedStiffness.ClearRows(index);
			simplifiedStiffness.ClearColumns(index);

			// Set the diagonal element to 1
			foreach (var i in index)
				simplifiedStiffness[i, i] = 1;

			if (!simplifyZeroRows)
				return simplifiedStiffness;

			// Verify rows
			for (var i = 0; i < simplifiedStiffness.RowCount; i++)
			{
				// Verify what line of the matrix is composed of zeroes
				if (simplifiedStiffness.Row(i).Exists(num => !num.ApproxZero()) && simplifiedStiffness.Column(i).Exists(num => !num.ApproxZero()))
					continue;

				// The row is composed of only zeroes, so the displacement must be zero
				// Set the diagonal element to 1
				simplifiedStiffness[i, i] = 1;
			}


			return simplifiedStiffness;
		}
	}
}