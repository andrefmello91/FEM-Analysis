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
	///     Stiffness matrix class.
	/// </summary>
	public class StiffnessMatrix : IUnitConvertible<ForcePerLengthUnit>, ICloneable<StiffnessMatrix>, IEquatable<StiffnessMatrix>
	{

		#region Fields

		private ForcePerLengthUnit _unit;

		/// <summary>
		///     The corresponding matrix, with components in <see cref="Unit" />.
		/// </summary>
		private Matrix<double> _value;

		#endregion

		#region Properties

		#region Interface Implementations

		/// <inheritdoc />
		public ForcePerLengthUnit Unit
		{
			get => _unit;
			set => ChangeUnit(value);
		}

		/// <summary>
		///		The index of constrained DoFs.
		/// </summary>
		public int[]? ConstraintIndex { get; set; }

		/// <inheritdoc cref="Matrix{T}.RowCount"/>
		public int Rows => _value.RowCount;
		
		/// <inheritdoc cref="Matrix{T}.ColumnCount"/>
		public int Columns => _value.ColumnCount;

		/// <summary>
		///		Get/set the value at these indexes.
		/// </summary>
		/// <param name="rowIndex">The row of the required element.</param>
		/// <param name="columnIndex">The column of the required element.</param>
		public ForcePerLength this[int rowIndex, int columnIndex]
		{
			get => ForcePerLength.From(_value[rowIndex, columnIndex], _unit);
			set => _value[rowIndex, columnIndex] = value.As(_unit);
		}

		#endregion

		#endregion

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
		{
			_value = value;
			_unit  = unit;
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public StiffnessMatrix Convert(ForcePerLengthUnit unit) => new(ForcePerLength.From(1, Unit).As(unit) * _value, unit);

		#region Interface Implementations

		/// <inheritdoc cref="Matrix{T}.Row(int)"/>
		public Vector<double> Row(int index) => _value.Row(index);
		
		/// <inheritdoc cref="Matrix{T}.Column(int)"/>
		public Vector<double> Column(int index) => _value.Column(index);

		/// <inheritdoc cref="Matrix{T}.ClearRows(int[])"/>
		public void ClearRows(params int[] indexes) => _value.ClearRows(indexes);
		
		/// <inheritdoc cref="Matrix{T}.ClearColumns(int[])"/>
		public void ClearColumns(params int[] indexes) => _value.ClearColumns(indexes);

		/// <inheritdoc cref="Matrix{T}.Clear"/>
		public void Clear() => _value.Clear();
		
		/// <inheritdoc cref="Matrix{T}.Determinant"/>
		public double Determinant() => _value.Determinant();
		
		/// <inheritdoc cref="Matrix{T}.Transpose()"/>
		public StiffnessMatrix Transpose() => new (_value.Transpose(), _unit);

		/// <summary>
		///		Get the simplified stiffness matrix by the constrained DoFs.
		/// </summary>
		public Matrix<double> Simplified() => ConstraintIndex is not null
			? SimplifiedStiffness(_value, ConstraintIndex)
			: _value;
		
		/// <summary>
		///		Transform this stiffness to another coordinate system.
		/// </summary>
		/// <param name="transformationMatrix">The transformation matrix.</param>
		/// <returns>
		///		A new <see cref="StiffnessMatrix"/> with transformed components.
		/// </returns>
		/// <exception cref="ArgumentException">If the dimensions of <paramref name="transformationMatrix"/> don't conform with this.</exception>
		public StiffnessMatrix Transform(Matrix<double> transformationMatrix)
		{
			var value = transformationMatrix.Transpose() * _value * transformationMatrix;
			
			return
				new StiffnessMatrix(value, _unit);
		}

		/// <inheritdoc />
		public void ChangeUnit(ForcePerLengthUnit unit)
		{
			if (_unit == unit)
				return;
			
			// Multiply matrix
			_value *= ForcePerLength.From(1, Unit).As(unit);
			
			// Set
			_unit = unit;
		}

		/// <inheritdoc />
		public StiffnessMatrix Clone() => new(_value.Clone(), _unit);

		/// <inheritdoc />
		public bool Equals(StiffnessMatrix? other) =>
			other is not null && _unit == other._unit && _value.Equals(other._value);

		/// <inheritdoc />
		IUnitConvertible<ForcePerLengthUnit> IUnitConvertible<ForcePerLengthUnit>.Convert(ForcePerLengthUnit unit) => Convert(unit);

		#endregion

		#region Object override

		/// <inheritdoc />
		public override bool Equals(object? obj) =>
			obj is StiffnessMatrix other && Equals(other);

		/// <inheritdoc />
		public override int GetHashCode() => (int) _unit * _value.GetHashCode();

		/// <inheritdoc />
		public override string ToString() =>
			$"Unit: {Unit} \n" +
			$"Value: {_value}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the corresponding <see cref="Matrix{T}" /> value of a <see cref="StiffnessMatrix" />.
		/// </summary>
		public static implicit operator Matrix<double>(StiffnessMatrix stiffnessMatrix) => stiffnessMatrix._value;

		/// <summary>
		///     Create a <see cref="StiffnessMatrix" /> from a <see cref="Matrix{T}" />, with components in
		///     <see cref="ForcePerLengthUnit.NewtonPerMillimeter" />.
		/// </summary>
		public static implicit operator StiffnessMatrix(Matrix<double> value) => new(value);

		/// <returns>
		///     True if objects are equal.
		/// </returns>
		public static bool operator ==(StiffnessMatrix left, StiffnessMatrix right) => left.IsEqualTo(right);

		/// <returns>
		///     True if objects are not equal.
		/// </returns>
		public static bool operator !=(StiffnessMatrix left, StiffnessMatrix right) => left.IsNotEqualTo(right);

		/// <returns>
		///     A new stiffness matrix with summed components in <paramref name="left" />'s unit.
		/// </returns>
		public static StiffnessMatrix operator +(StiffnessMatrix left, StiffnessMatrix right) => new(left._value + right.Convert(left.Unit)._value, left.Unit);

		/// <returns>
		///     A new stiffness matrix with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		public static StiffnessMatrix operator -(StiffnessMatrix left, StiffnessMatrix right) => new(left._value - right.Convert(left.Unit)._value, left.Unit);

		/// <returns>
		///     A new stiffness matrix with components multiplied by a value
		/// </returns>
		public static StiffnessMatrix operator *(double value, StiffnessMatrix right) => new(value * right._value, right.Unit);

		/// <inheritdoc cref="op_Multiply(double, StiffnessMatrix) " />
		public static StiffnessMatrix operator *(StiffnessMatrix left, double value) => value * left;

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