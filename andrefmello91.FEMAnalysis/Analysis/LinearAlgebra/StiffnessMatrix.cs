using System;
using System.Collections;
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
		public List<int>? ConstraintIndex { get; set; }

		/// <inheritdoc cref="Matrix{T}.RowCount"/>
		public int Rows => Values.GetLength(0);
		
		/// <inheritdoc cref="Matrix{T}.ColumnCount"/>
		public int Columns => Values.GetLength(1);

		/// <summary>
		///		Get/set the value at these indexes.
		/// </summary>
		/// <param name="rowIndex">The row of the required element.</param>
		/// <param name="columnIndex">The column of the required element.</param>
		public TQuantity this[int rowIndex, int columnIndex]
		{
			get => (TQuantity) Values[rowIndex, columnIndex].As(_unit);
			set => Values[rowIndex, columnIndex] = value.As(_unit);
		}

		/// <summary>
		///     The corresponding matrix, with components in <see cref="Unit" />.
		/// </summary>
		protected readonly double[,] Values;

		#endregion

		#endregion

		#region Constructors

		/// <summary>
		///     Create a stiffness matrix.
		/// </summary>
		/// <param name="values">The array of values.</param>
		/// <param name="unit">The unit of <paramref name="values" />'s components</param>
		public StiffnessMatrix(double[,] values, TUnit unit)
		{
			Values = (double[,]) values.Clone();
			_unit = unit;
		}
		
		/// <inheritdoc cref="StiffnessMatrix{T,T}(double[,], TUnit)" />
		public StiffnessMatrix(Matrix<double> value, TUnit unit)
		{
			Values = value.ToArray();
			_unit  = unit;
		}

		/// <inheritdoc cref="StiffnessMatrix{T,T}(double[,], TUnit)" />
		public StiffnessMatrix(TQuantity[,] value)
		{
			Values = value.GetValues<TQuantity, TUnit>();
			_unit  = value[0, 0].Unit;
		}

		#endregion

		#region Methods

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public StiffnessMatrix<TQuantity, TUnit> Convert(TUnit unit) => Unit.Equals(unit)
			? Clone()
			: new StiffnessMatrix<TQuantity, TUnit>(Values, unit)
			{
				ConstraintIndex = ConstraintIndex
			};

		#region Interface Implementations

		/// <inheritdoc cref="Matrix{T}.Row(int)"/>
		public Vector<double> Row(int index) => Values
			.GetRow(index)
			.ToVector();

		/// <inheritdoc cref="Matrix{T}.Column(int)"/>
		public Vector<double> Column(int index) => Values
			.GetColumn(index)
			.ToVector();

		/// <inheritdoc cref="Matrix{T}.ClearRows(int[])"/>
		public void ClearRows(params int[] indexes)
		{
			foreach (var i in indexes)
				for (var j = 0; j < Columns; j++)
					Values[i, j] = 0;
		}

		/// <inheritdoc cref="Matrix{T}.ClearColumns(int[])"/>
		public void ClearColumns(params int[] indexes)
		{
			foreach (var j in indexes) 
				for (var i = 0; i < Rows; i++)
					Values[i, j] = 0;
		}

		/// <inheritdoc cref="Matrix{T}.Clear"/>
		public void Clear()
		{
			for (var i = 0; i < Rows; i++)
			for (var j = 0; j < Columns; j++)
				Values[i, j] = 0;
		}
		
		/// <inheritdoc cref="Matrix{T}.Determinant"/>
		public double Determinant() => Values
			.ToMatrix()
			.Determinant();
		
		/// <inheritdoc cref="Matrix{T}.Transpose()"/>
		public StiffnessMatrix<TQuantity, TUnit> Transpose() => new (Values.ToMatrix().Transpose(), Unit);

		/// <summary>
		///		Get the simplified stiffness matrix by the constrained DoFs.
		/// </summary>
		/// <param name="threshold">A value for setting all values whose absolute value is smaller than to zero. If null, this is not applied.</param>
		/// <returns>
		///		The simplified <see cref="Matrix{T}"/>.
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

		/// <inheritdoc cref="Simplified(double?)"/>
		public Matrix<double> Simplified(TQuantity? threshold) => Simplified(threshold?.As(Unit));
		
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
			var value = transformationMatrix.Transpose() * Values.ToMatrix() * transformationMatrix;
			
			return
				new StiffnessMatrix<TQuantity, TUnit>(value, Unit)
				{
					ConstraintIndex = ConstraintIndex
				};
		}

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
		public StiffnessMatrix<TQuantity, TUnit> Clone() => new (Values, Unit)
		{
			ConstraintIndex = ConstraintIndex
		};

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

		/// <returns>
		///     A new stiffness matrix with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentOutOfRangeException">If left and right don't have the same dimensions.</exception>
		public static StiffnessMatrix<TQuantity, TUnit> operator +(StiffnessMatrix<TQuantity, TUnit> left, StiffnessMatrix<TQuantity, TUnit> right) =>
			new(left.Values.ToMatrix() + right.Convert(left.Unit).Values.ToMatrix(), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};

		/// <returns>
		///     A new stiffness matrix with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentOutOfRangeException">If left and right don't have the same dimensions.</exception>
		public static StiffnessMatrix<TQuantity, TUnit> operator -(StiffnessMatrix<TQuantity, TUnit> left, StiffnessMatrix<TQuantity, TUnit> right) =>
			new(left.Values.ToMatrix() - right.Convert(left.Unit).Values.ToMatrix(), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};


		/// <returns>
		///     A new stiffness matrix with components multiplied by a value
		/// </returns>
		public static StiffnessMatrix<TQuantity, TUnit> operator *(double multiplier, StiffnessMatrix<TQuantity, TUnit> matrix) => new(multiplier * matrix.Values.ToMatrix(), matrix.Unit)
		{
			ConstraintIndex = matrix.ConstraintIndex
		};


		/// <inheritdoc cref="op_Multiply(double, StiffnessMatrix{TQuantity,TUnit}) " />
		public static StiffnessMatrix<TQuantity, TUnit> operator *(StiffnessMatrix<TQuantity, TUnit> matrix, double multiplier) => multiplier * matrix;
		
		/// <inheritdoc cref="Matrix{T}.op_Division(Matrix{T}, T)"/>
		public static StiffnessMatrix<TQuantity, TUnit> operator / (StiffnessMatrix<TQuantity, TUnit> matrix, double divisor) => new(matrix.Values.ToMatrix() / divisor, matrix.Unit)
		{
			ConstraintIndex = matrix.ConstraintIndex
		};


		/// <inheritdoc cref="Matrix{T}.op_UnaryNegation"/>
		public static StiffnessMatrix<TQuantity, TUnit> operator -(StiffnessMatrix<TQuantity, TUnit> matrix) => new (-matrix.Values.ToMatrix(), matrix.Unit)
		{
			ConstraintIndex = matrix.ConstraintIndex
		};


		#endregion

	}
	
	/// <summary>
	///     Default stiffness matrix class.
	/// </summary>
	/// <remarks>
	///		Unit is <see cref="ForcePerLengthUnit"/>.
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

		///  <summary>
		/// 		Solve a system d = K f.
		///  </summary>
		///  <param name="forceVector"></param>
		///  <param name="useSimplified">Use simplified matrix and vector?</param>
		///  <returns>
		/// 		The resulting displacement vector, with components in <see cref="LengthUnit.Millimeter"/>.
		///  </returns>
		public DisplacementVector Solve(ForceVector forceVector, bool useSimplified = true)
		{
			// Convert
			var stiffnessMatrix = Unit is ForcePerLengthUnit.NewtonPerMillimeter
				? this
				: Convert(ForcePerLengthUnit.NewtonPerMillimeter);

			var forces = forceVector.Unit is ForceUnit.Newton
				? forceVector
				: forceVector.Convert(ForceUnit.Newton);

			Matrix<double> k = useSimplified
				? stiffnessMatrix.Simplified()
				: stiffnessMatrix;

			var f = useSimplified
				? forces
				: forces.Simplified();

			// Solve
			var d = k.Solve(f);

			return
				new DisplacementVector(d);
		}
		
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
		///		Create a stiffness matrix with zero elements.
		/// </summary>
		/// <param name="size">The size of the matrix.</param>
		public static StiffnessMatrix Zero(int size) => new (Matrix<double>.Build.Dense(size, size));
		
		/// <summary>
		///     Assemble the global stiffness matrix.
		/// </summary>
		/// <param name="femInput">The <see cref="FEMInput" /></param>
		public static StiffnessMatrix Assemble(IFEMInput femInput)
		{
			var stiffness = Zero(femInput.NumberOfDoFs);
			stiffness.ConstraintIndex = femInput.ConstraintIndex;
			
			foreach (var element in femInput)
			{
				var dofIndex = element.DoFIndex;

				for (var i = 0; i < dofIndex.Length; i++)
				{
					// Global index
					var k = dofIndex[i];

					for (var j = 0; j < dofIndex.Length; j++)
					{
						// Global index
						var l = dofIndex[j];

						stiffness[k, l] += element.Stiffness[i, j];
					}
				}
			}

			return stiffness;
		}

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
		
		/// <summary>
		///		Create a force vector by multiplying the stiffness and the displacement vector.
		/// </summary>
		/// <remarks>
		///		This uses the simplified stiffness matrix.
		/// </remarks>
		/// <returns>
		///		The <see cref="ForceVector"/> with components in <see cref="ForceUnit.Newton"/>.
		/// </returns>
		public static ForceVector operator *(StiffnessMatrix stiffnessMatrix, DisplacementVector displacementVector)
		{
			// Convert
			var k = stiffnessMatrix.Convert(ForcePerLengthUnit.NewtonPerMillimeter).Simplified();
			var d = displacementVector.Convert(LengthUnit.Millimeter);
			
			// Multiply
			var f = k * d;

			return
				new ForceVector(f);
		}

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.op_Addition" />
		public static StiffnessMatrix operator +(StiffnessMatrix left, StiffnessMatrix right) => (StiffnessMatrix) ((StiffnessMatrix<ForcePerLength, ForcePerLengthUnit>) left + right);

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.op_Subtraction" />
		public static StiffnessMatrix operator -(StiffnessMatrix left, StiffnessMatrix right) => (StiffnessMatrix) ((StiffnessMatrix<ForcePerLength, ForcePerLengthUnit>) left - right);

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.op_Multiply(double, StiffnessMatrix{TQuantity,TUnit}) " />
		public static StiffnessMatrix operator *(double value, StiffnessMatrix right) => (StiffnessMatrix) (value * (StiffnessMatrix<ForcePerLength, ForcePerLengthUnit>) right);

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.op_Multiply(double, StiffnessMatrix{TQuantity,TUnit}) " />
		public static StiffnessMatrix operator *(StiffnessMatrix left, double value) => value * left;
		
		/// <inheritdoc cref="Matrix{T}.op_UnaryNegation"/>
		public static StiffnessMatrix operator -(StiffnessMatrix right) => new (-right.Values.ToMatrix(), right.Unit)
		{
			ConstraintIndex = right.ConstraintIndex
		};


		/// <inheritdoc cref="Matrix{T}.op_Division(Matrix{T}, T)"/>
		public static StiffnessMatrix operator / (StiffnessMatrix left, double value) => new(left.Values.ToMatrix() / value, left.Unit)
		{
			ConstraintIndex = left.ConstraintIndex
		};
		
		///  <summary>
		///      Calculate the stiffness increment for nonlinear analysis.
		///  </summary>
		///  <param name="solver">The nonlinear solver.</param>
		///  <returns>
		/// 		The <see cref="StiffnessMatrix"/> increment with current unit.
		///  </returns>
		/// <inheritdoc cref="TangentIncrement"/>
		public static StiffnessMatrix StiffnessIncrement(IIteration currentIteration, IIteration lastIteration, NonLinearSolver solver = NonLinearSolver.NewtonRaphson) =>
			solver switch
			{
				NonLinearSolver.Secant => SecantIncrement(currentIteration, lastIteration),
				_                      => TangentIncrement(currentIteration, lastIteration)
			};
		
		/// <summary>
		///     Calculate the secant stiffness increment.
		/// </summary>
		/// <inheritdoc cref="TangentIncrement(andrefmello91.FEMAnalysis.IIteration, andrefmello91.FEMAnalysis.IIteration)"/>
		private static StiffnessMatrix SecantIncrement(IIteration currentIteration, IIteration lastIteration)
		{
			// Calculate the variation of displacements and residual as vectors
			Vector<double>
				dU = (currentIteration.Displacements - lastIteration.Displacements).Convert(LengthUnit.Millimeter),
				dR = (currentIteration.ResidualForces - lastIteration.ResidualForces).Convert(ForceUnit.Newton);

			Matrix<double> stiffness = lastIteration.Stiffness.Convert(ForcePerLengthUnit.NewtonPerMillimeter);
			
			var inc = ((dR - stiffness * dU) / dU.Norm(2)).ToColumnMatrix() * dU.ToRowMatrix();

			return
				new StiffnessMatrix(inc)
				{
					ConstraintIndex = currentIteration.Stiffness.ConstraintIndex
				}
					.Convert(currentIteration.Stiffness.Unit);
		}

		/// <summary>
		///     Calculate the tangent stiffness increment.
		/// </summary>
		/// <param name="currentIteration">The current iteration.</param>
		/// <param name="lastIteration">The last solved iteration.</param>
		/// <returns>
		///		The <see cref="andrefmello91.FEMAnalysis.StiffnessMatrix"/> increment with current unit.
		/// </returns>
		private static StiffnessMatrix TangentIncrement(IIteration currentIteration, IIteration lastIteration)
		{
			// Get variations
			var dF = (currentIteration.InternalForces - lastIteration.InternalForces).Convert(ForceUnit.Newton);
			var dU = (currentIteration.Displacements - lastIteration.Displacements).Convert(LengthUnit.Millimeter);

			var inc = dF.ToColumnMatrix() * dU.ToRowMatrix();
			
			return
				new StiffnessMatrix(inc)
				{
					ConstraintIndex = currentIteration.Stiffness.ConstraintIndex
				}
					.Convert(currentIteration.Stiffness.Unit);
		}
	}
}