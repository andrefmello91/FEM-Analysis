using System;
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
	///     Generic stiffness matrix class.
	/// </summary>
	/// <typeparam name="TQuantity">The quantity that represents the value of components of the matrix.</typeparam>
	/// <typeparam name="TUnit">The unit enumeration that represents the quantity of the components of the matrix.</typeparam>
	public abstract class QuantityMatrix<TQuantity, TUnit> : DenseMatrix, IUnitConvertible<TUnit>, ICloneable<QuantityMatrix<TQuantity, TUnit>>, IEquatable<QuantityMatrix<TQuantity, TUnit>> where TQuantity : IQuantity<TUnit>
		where TUnit : Enum
	{

		#region Fields

		private TUnit _unit;

		#endregion

		#region Properties

		/// <summary>
		///     The index of constrained DoFs.
		/// </summary>
		public List<int>? ConstraintIndex { get; set; }

		/// <summary>
		///     Get/set the value at these indexes.
		/// </summary>
		/// <param name="rowIndex">The row of the required element.</param>
		/// <param name="columnIndex">The column of the required element.</param>
		public new TQuantity this[int rowIndex, int columnIndex]
		{
			get => (TQuantity) base[rowIndex, columnIndex].As(_unit);
			set => base[rowIndex, columnIndex] = value.As(_unit);
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
		///     Create a quantity matrix.
		/// </summary>
		/// <inheritdoc />
		/// <param name="unit">The unit of matrix components</param>
		protected QuantityMatrix(DenseColumnMajorMatrixStorage<double> storage, TUnit unit)
			: base(storage) =>
			_unit = unit;

		/// <summary>
		///     Create a quantity matrix.
		/// </summary>
		/// <inheritdoc />
		/// <param name="unit">The unit of <paramref name="values" />'s components</param>
		protected QuantityMatrix(int rows, int columns, TUnit unit)
			: base(rows, columns) =>
			_unit = unit;

		/// <summary>
		///     Create a quantity matrix.
		/// </summary>
		/// <inheritdoc />
		/// <param name="unit">The unit of <paramref name="values" />'s components</param>
		protected QuantityMatrix(int rows, int columns, double[] storage, TUnit unit)
			: base(rows, columns, storage) =>
			_unit = unit;

		/// <summary>
		///     Create a quantity matrix.
		/// </summary>
		/// <param name="values">The array of values.</param>
		/// <param name="unit">The unit of <paramref name="values" />'s components</param>
		protected QuantityMatrix(double[,] values, TUnit unit)
			: this(DenseColumnMajorMatrixStorage<double>.OfArray(values), unit)
		{
		}

		/// <inheritdoc cref="QuantityMatrix{TQuantity,TUnit}(double[,],TUnit)" />
		protected QuantityMatrix(Matrix<double> value, TUnit unit)
			: this(value.ToArray(), unit)
		{
		}

		/// <inheritdoc cref="QuantityMatrix{TQuantity,TUnit}(double[,],TUnit)" />
		protected QuantityMatrix(TQuantity[,] value)
			: this(value.GetValues<TQuantity, TUnit>(), value[0, 0].Unit)
		{
		}

		#endregion

		#region Methods

		/// <summary>
		///     Get the global stiffness simplified.
		/// </summary>
		/// <param name="stiffness">The global stiffness <see cref="Matrix{T}" /> to simplify.</param>
		/// <param name="indexes">The DoF indexes to simplify matrix.</param>
		/// <param name="simplifyZeroRows">Simplify matrix at rows containing only zero elements?</param>
		protected static Matrix<double> SimplifiedStiffness(Matrix<double> stiffness, IEnumerable<int> indexes, bool simplifyZeroRows = true)
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

		/// <inheritdoc cref="Matrix{T}.Add(T)" />
		public QuantityMatrix<TQuantity, TUnit> Add(TQuantity scalar)
		{
			var result = Clone();

			if (scalar.Value.Equals(0))
				return result;

			result.Clear();

			DoAdd(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Matrix{T}.Add(Matrix{T})" />
		public QuantityMatrix<TQuantity, TUnit> Add(QuantityMatrix<TQuantity, TUnit> other)
		{
			var result = other.Clone();

			result.Clear();

			Add(other, result);

			return result;
		}

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public abstract QuantityMatrix<TQuantity, TUnit> Convert(TUnit unit);

		/// <inheritdoc cref="Matrix{T}.Divide(T)" />
		public new QuantityMatrix<TQuantity, TUnit> Divide(double scalar)
		{
			var result = Clone();

			if (scalar.Equals(1))
				return result;

			result.Clear();

			DoDivide(scalar, result);

			return result;
		}

		/// <inheritdoc cref="Matrix{T}.DivideByThis(T)" />
		public Matrix<double> DivideByThis(TQuantity scalar)
		{
			var result = Build.Dense(RowCount, ColumnCount);

			if (scalar.Value.Equals(0))
				return result;

			DoDivideByThis(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Matrix{T}.Multiply(T)" />
		public new QuantityMatrix<TQuantity, TUnit> Multiply(double scalar)
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

		/// <inheritdoc cref="Matrix{T}.Multiply(Matrix{T})" />
		public new QuantityMatrix<TQuantity, TUnit> Multiply(Matrix<double> other)
		{
			var result = Clone();

			result.Clear();

			Multiply(other, result);

			return result;
		}

		/// <summary>
		///     Multiply the other matrix by this.
		/// </summary>
		/// <param name="other">The matrix to multiply by this</param>
		public QuantityMatrix<TQuantity, TUnit> MultiplyBy(Matrix<double> other)
		{
			var result = Clone();

			result.Clear();

			other.Multiply(this, result);

			return result;
		}

		/// <inheritdoc cref="Matrix{T}.Negate()" />
		public new QuantityMatrix<TQuantity, TUnit> Negate()
		{
			var result = Clone();

			result.Clear();

			DoNegate(result);

			return result;
		}

		/// <summary>
		///     Get the simplified stiffness matrix by the constrained DoFs.
		/// </summary>
		/// <remarks>
		///     This uses the default tolerance.
		/// </remarks>
		/// <returns>
		///     The simplified <see cref="Matrix{T}" />.
		/// </returns>
		public abstract QuantityMatrix<TQuantity, TUnit> Simplified();

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
		public abstract QuantityMatrix<TQuantity, TUnit> Simplified(double? threshold);

		/// <inheritdoc cref="Simplified(double?)" />
		public QuantityMatrix<TQuantity, TUnit> Simplified(TQuantity? threshold) => Simplified(threshold?.As(Unit));

		/// <inheritdoc cref="Matrix{T}.Subtract(T)" />
		public QuantityMatrix<TQuantity, TUnit> Subtract(TQuantity scalar)
		{
			var result = Clone();

			if (scalar.Value.Equals(0))
				return result;

			result.Clear();

			DoSubtract(scalar.As(Unit), result);

			return result;
		}

		/// <inheritdoc cref="Matrix{T}.Subtract(Matrix{T})" />
		public QuantityMatrix<TQuantity, TUnit> Subtract(QuantityMatrix<TQuantity, TUnit> other)
		{
			var result = other.Clone();

			result.Clear();

			Subtract(other, result);

			return result;
		}

		/// <inheritdoc cref="Matrix{T}.SubtractFrom(T)" />
		public QuantityMatrix<TQuantity, TUnit> SubtractFrom(TQuantity scalar)
		{
			var result = Clone();

			if (scalar.Value.Equals(0))
				return result;

			result.Clear();

			DoSubtractFrom(scalar.As(Unit), result);

			return result;
		}

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
		public abstract QuantityMatrix<TQuantity, TUnit> Transform(Matrix<double> transformationMatrix);

		/// <inheritdoc cref="Matrix{T}.Transpose()" />
		public abstract QuantityMatrix<TQuantity, TUnit> Transpose();

		#region Interface Implementations

		/// <inheritdoc />
		public abstract QuantityMatrix<TQuantity, TUnit> Clone();

		/// <inheritdoc />
		public bool Equals(QuantityMatrix<TQuantity, TUnit>? other) =>
			other is not null &&
			Equals(other.Unit.Equals(Unit)
				? other
				: other.Convert(Unit));

		/// <inheritdoc />
		public void ChangeUnit(TUnit unit)
		{
			if (_unit.Equals(unit))
				return;

			// Multiply matrix
			for (var i = 0; i < RowCount; i++)
			for (var j = 0; j < ColumnCount; j++)
				base[i, j] = this[i, j].As(unit);

			// Set
			_unit = unit;
		}

		/// <inheritdoc />
		IUnitConvertible<TUnit> IUnitConvertible<TUnit>.Convert(TUnit unit) => Convert(unit);

		#endregion

		#region Object override

		/// <inheritdoc />
		public override bool Equals(object? obj) =>
			obj is QuantityMatrix<TQuantity, TUnit> other && Equals(other);

		/// <inheritdoc />
		public override int GetHashCode() => _unit.GetHashCode() * Values.GetHashCode();

		/// <inheritdoc cref="Add(TQuantity)" />
		public static QuantityMatrix<TQuantity, TUnit> operator +(QuantityMatrix<TQuantity, TUnit> left, TQuantity right) => left.Add(right);

		/// <inheritdoc cref="Add(TQuantity)" />
		public static QuantityMatrix<TQuantity, TUnit> operator +(TQuantity left, QuantityMatrix<TQuantity, TUnit> right) => right.Add(left);

		/// <inheritdoc cref="Add(QuantityMatrix{TQuantity, TUnit})" />
		public static QuantityMatrix<TQuantity, TUnit> operator +(QuantityMatrix<TQuantity, TUnit> left, QuantityMatrix<TQuantity, TUnit> right) => left.Add(right);

		/// <inheritdoc cref="Divide(double)" />
		public static QuantityMatrix<TQuantity, TUnit> operator /(QuantityMatrix<TQuantity, TUnit> matrix, double divisor) => matrix.Divide(divisor);

		/// <inheritdoc cref="DivideByThis(TQuantity)" />
		public static Matrix<double> operator /(TQuantity quantity, QuantityMatrix<TQuantity, TUnit> matrix) => matrix.DivideByThis(quantity);

		/// <returns>
		///     True if objects are equal.
		/// </returns>
		public static bool operator ==(QuantityMatrix<TQuantity, TUnit>? left, QuantityMatrix<TQuantity, TUnit>? right) => left.IsEqualTo(right);

		/// <returns>
		///     True if objects are not equal.
		/// </returns>
		public static bool operator !=(QuantityMatrix<TQuantity, TUnit>? left, QuantityMatrix<TQuantity, TUnit>? right) => left.IsNotEqualTo(right);

		/// <inheritdoc cref="Multiply(double)" />
		public static QuantityMatrix<TQuantity, TUnit> operator *(QuantityMatrix<TQuantity, TUnit> matrix, double multiplier) => matrix.Multiply(multiplier);

		/// <inheritdoc cref="Multiply(double)" />
		public static QuantityMatrix<TQuantity, TUnit> operator *(double multiplier, QuantityMatrix<TQuantity, TUnit> matrix) => matrix.Multiply(multiplier);

		/// <inheritdoc cref="Multiply(Matrix{double})" />
		public static QuantityMatrix<TQuantity, TUnit> operator *(QuantityMatrix<TQuantity, TUnit> matrix, Matrix<double> multiplier) => matrix.Multiply(multiplier);

		/// <inheritdoc cref="MultiplyBy(Matrix{double})" />
		public static QuantityMatrix<TQuantity, TUnit> operator *(Matrix<double> multiplier, QuantityMatrix<TQuantity, TUnit> matrix) => matrix.MultiplyBy(multiplier);

		/// <inheritdoc cref="Subtract(TQuantity)" />
		public static QuantityMatrix<TQuantity, TUnit> operator -(QuantityMatrix<TQuantity, TUnit> left, TQuantity right) => left.Subtract(right);

		/// <inheritdoc cref="Subtract(QuantityMatrix{TQuantity, TUnit})" />
		public static QuantityMatrix<TQuantity, TUnit> operator -(QuantityMatrix<TQuantity, TUnit> left, QuantityMatrix<TQuantity, TUnit> right) => left.Subtract(right);

		/// <inheritdoc cref="SubtractFrom(TQuantity)" />
		public static QuantityMatrix<TQuantity, TUnit> operator -(TQuantity left, QuantityMatrix<TQuantity, TUnit> right) => right.SubtractFrom(left);

		/// <inheritdoc cref="Negate" />
		public static QuantityMatrix<TQuantity, TUnit> operator -(QuantityMatrix<TQuantity, TUnit> right) => right.Negate();

		/// <inheritdoc cref="object.ToString" />
		public new string ToString() =>
			$"Unit: {Unit} \n" +
			$"Value: {base.ToString()}";

		#endregion

		#endregion

	}
}