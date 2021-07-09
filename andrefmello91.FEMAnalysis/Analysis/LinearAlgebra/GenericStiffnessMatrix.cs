using System;
using System.Collections.Generic;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Generic stiffness matrix class.
	/// </summary>
	/// <typeparam name="TQuantity">The quantity that represents the value of components of the matrix.</typeparam>
	/// <typeparam name="TUnit">The unit enumeration that represents the quantity of the components of the matrix.</typeparam>
	public abstract class StiffnessMatrix<TQuantity, TUnit> : IUnitConvertible<TUnit>, ICloneable<StiffnessMatrix<TQuantity, TUnit>>, IEquatable<StiffnessMatrix<TQuantity, TUnit>>
		where TQuantity : IQuantity<TUnit>
		where TUnit : Enum
	{

		#region Fields

		/// <summary>
		///     The corresponding matrix, with components in <see cref="Unit" />.
		/// </summary>
		protected readonly double[,] Values;

		private TUnit _unit;

		#endregion

		#region Properties

		/// <inheritdoc cref="Matrix{T}.ColumnCount" />
		public int Columns => Values.GetLength(1);

		/// <summary>
		///     The index of constrained DoFs.
		/// </summary>
		public List<int>? ConstraintIndex { get; set; }

		/// <summary>
		///     Get/set the value at these indexes.
		/// </summary>
		/// <param name="rowIndex">The row of the required element.</param>
		/// <param name="columnIndex">The column of the required element.</param>
		public TQuantity this[int rowIndex, int columnIndex]
		{
			get => (TQuantity) Values[rowIndex, columnIndex].As(_unit);
			set => Values[rowIndex, columnIndex] = value.As(_unit);
		}

		/// <inheritdoc cref="Matrix{T}.RowCount" />
		public int Rows => Values.GetLength(0);

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
		///     Create a stiffness matrix.
		/// </summary>
		/// <param name="values">The array of values.</param>
		/// <param name="unit">The unit of <paramref name="values" />'s components</param>
		protected StiffnessMatrix(double[,] values, TUnit unit)
		{
			Values = (double[,]) values.Clone();
			_unit  = unit;
		}

		/// <inheritdoc cref="StiffnessMatrix{T,T}(double[,], TUnit)" />
		protected StiffnessMatrix(Matrix<double> value, TUnit unit)
		{
			Values = value.ToArray();
			_unit  = unit;
		}

		/// <inheritdoc cref="StiffnessMatrix{T,T}(double[,], TUnit)" />
		protected StiffnessMatrix(TQuantity[,] value)
		{
			Values = value.GetValues<TQuantity, TUnit>();
			_unit  = value[0, 0].Unit;
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="Matrix{T}.Clear" />
		public void Clear()
		{
			for (var i = 0; i < Rows; i++)
			for (var j = 0; j < Columns; j++)
				Values[i, j] = 0;
		}

		/// <inheritdoc cref="Matrix{T}.ClearColumns(int[])" />
		public void ClearColumns(params int[] indexes)
		{
			foreach (var j in indexes)
				for (var i = 0; i < Rows; i++)
					Values[i, j] = 0;
		}

		/// <inheritdoc cref="Matrix{T}.ClearRows(int[])" />
		public void ClearRows(params int[] indexes)
		{
			foreach (var i in indexes)
				for (var j = 0; j < Columns; j++)
					Values[i, j] = 0;
		}

		/// <inheritdoc cref="Matrix{T}.Column(int)" />
		public Vector<double> Column(int index) => Values
			.GetColumn(index)
			.ToVector();

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public abstract StiffnessMatrix<TQuantity, TUnit> Convert(TUnit unit);

		/// <inheritdoc cref="Matrix{T}.Determinant" />
		public double Determinant() => Values
			.ToMatrix()
			.Determinant();

		/// <inheritdoc cref="Matrix{T}.Row(int)" />
		public Vector<double> Row(int index) => Values
			.GetRow(index)
			.ToVector();

		/// <summary>
		///     Get the simplified stiffness matrix by the constrained DoFs.
		/// </summary>
		/// <param name="threshold">
		///     A value for setting all values whose absolute value is smaller than to zero. If null, this is
		///     not applied.
		/// </param>
		/// <returns>
		///     The simplified <see cref="Matrix{T}" />.
		/// </returns>
		public Matrix<double> Simplified(double? threshold = null)
		{
			var value = Values.ToMatrix();

			if (threshold.HasValue)
				value.CoerceZero(threshold.Value);

			return ConstraintIndex is not null
				? StiffnessMatrix.SimplifiedStiffness(value, ConstraintIndex)
				: value;
		}

		/// <inheritdoc cref="Simplified(double?)" />
		public Matrix<double> Simplified(TQuantity? threshold) => Simplified(threshold?.As(Unit));

		/// <summary>
		///     Transform this stiffness to another coordinate system.
		/// </summary>
		/// <param name="transformationMatrix">The transformation matrix.</param>
		/// <returns>
		///     A new <see cref="StiffnessMatrix" /> with transformed components.
		/// </returns>
		/// <exception cref="ArgumentException">
		///     If the dimensions of <paramref name="transformationMatrix" /> don't conform with
		///     this.
		/// </exception>
		public abstract StiffnessMatrix<TQuantity, TUnit> Transform(Matrix<double> transformationMatrix);

		/// <inheritdoc cref="Matrix{T}.Transpose()" />
		public abstract StiffnessMatrix<TQuantity, TUnit> Transpose();

		#region Interface Implementations

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;

			// Multiply matrix
			for (var i = 0; i < Rows; i++)
			for (var j = 0; j < Columns; j++)
				Values[i, j] = this[i, j].As(unit);

			// Set
			_unit = unit;
		}

		/// <inheritdoc />
		public abstract StiffnessMatrix<TQuantity, TUnit> Clone();

		/// <inheritdoc />
		public bool Equals(StiffnessMatrix<TQuantity, TUnit>? other) =>
			other is not null &&
			Values.ToMatrix().Equals(other.Convert(Unit).Values.ToMatrix());

		/// <inheritdoc />
		IUnitConvertible<TUnit> IUnitConvertible<TUnit>.Convert(TUnit unit) => Convert(unit);

		#endregion

		#region Object override

		/// <inheritdoc />
		public override bool Equals(object? obj) =>
			obj is StiffnessMatrix<TQuantity, TUnit> other && Equals(other);

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
		///     Get the corresponding <see cref="Matrix{T}" /> value of a <see cref="StiffnessMatrix" />.
		/// </summary>
		public static implicit operator Matrix<double>(StiffnessMatrix<TQuantity, TUnit> stiffnessMatrix) => stiffnessMatrix.Values.ToMatrix();

		/// <returns>
		///     True if objects are equal.
		/// </returns>
		public static bool operator ==(StiffnessMatrix<TQuantity, TUnit>? left, StiffnessMatrix<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <returns>
		///     True if objects are not equal.
		/// </returns>
		public static bool operator !=(StiffnessMatrix<TQuantity, TUnit>? left, StiffnessMatrix<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		#endregion

	}
}