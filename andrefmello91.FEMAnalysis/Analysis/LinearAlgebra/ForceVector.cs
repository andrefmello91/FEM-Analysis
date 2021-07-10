using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Force vector class.
	/// </summary>
	/// <remarks>
	///     Unit is <see cref="ForceUnit" />.
	///     <para>
	///         Quantity is <see cref="Force" />.
	///     </para>
	/// </remarks>
	public class ForceVector : QuantityVector<Force, ForceUnit>
	{

		/// <summary>
		///     Default tolerance for force vector.
		/// </summary>
		private static Force Tolerance { get; } = PlaneForce.Tolerance;

		#region Constructors

		/// <inheritdoc />
		/// <remarks>
		///     Default unit is Newton.
		/// </remarks>
		public ForceVector(IEnumerable<double> values, ForceUnit unit = ForceUnit.Newton)
			: base(values, unit)
		{
		}

		/// <inheritdoc />
		public ForceVector(IEnumerable<Force> values)
			: base(values)
		{
		}

		#endregion

		#region Methods

		/// <summary>
		///     Assemble the global external force vector.
		/// </summary>
		/// <param name="femInput">Finite element input.</param>
		public static ForceVector AssembleExternal(IFEMInput femInput)
		{
			// Initialize the force vector
			var f = Zero(femInput.NumberOfDoFs);

			// Read the nodes data
			foreach (var grip in femInput.Grips)
			{
				// Get DoF indexes
				var index = grip.DoFIndex;
				int
					i = index[0],
					j = index[1];

				// Set to force vector
				f[i] = grip.Force.X;
				f[j] = grip.Force.Y;
			}

			return f;
		}

		/// <summary>
		///     Assemble the global internal force vector.
		/// </summary>
		/// <inheritdoc cref="AssembleExternal" />
		public static ForceVector AssembleInternal(IFEMInput femInput)
		{
			var fi = Zero(femInput.NumberOfDoFs);

			foreach (var element in femInput)
			{
				var dofIndex = element.DoFIndex;

				for (var i = 0; i < dofIndex.Length; i++)
				{
					// DoF index
					var j = dofIndex[i];

					// Add values
					fi[j] += element.Forces[i];
				}
			}

			return fi;
		}

		/// <summary>
		///     Create a force vector with zero elements.
		/// </summary>
		/// <param name="size">The size of the vector.</param>
		public new static ForceVector Zero(int size) => new(new double[size]);


		/// <inheritdoc />
		public override QuantityVector<Force, ForceUnit> Convert(ForceUnit unit) =>Unit == unit 
			? Clone() 
			: new ForceVector(Values.GetQuantities<Force, ForceUnit>(Unit).GetValues(unit), unit);

		/// <inheritdoc />
		public override QuantityVector<Force, ForceUnit> Simplified(IEnumerable<int>? indexes = null) => Simplified(Tolerance, indexes);

		/// <inheritdoc />
		public override QuantityVector<Force, ForceUnit> Simplified(double? threshold, IEnumerable<int>? indexes = null)
		{
			var simplified = Values.ToVector();

			if (indexes is not null)
				foreach (var index in indexes)
					simplified[index] = 0;

			if (threshold.HasValue)
				simplified.CoerceZero(threshold.Value);

			return
				new ForceVector(simplified, Unit);
		}

		/// <inheritdoc cref="ICloneable.Clone" />
		public override QuantityVector<Force, ForceUnit> Clone() => new ForceVector(Values, Unit);


		#endregion

		#region Operators

		/// <returns>
		///     A new vector with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ForceVector operator +(ForceVector left, ForceVector right)
		{
			right = right.Unit == left.Unit
				? right
				: (ForceVector) right.Convert(left.Unit);

			var vec = (Vector<double>) left + right;

			return
				new ForceVector(vec, left.Unit);
		}

		/// <returns>
		///     A new vector with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ForceVector operator -(ForceVector left, ForceVector right) 
		{
			right = right.Unit == left.Unit
				? right
				: (ForceVector) right.Convert(left.Unit);

			var vec = (Vector<double>) left - right;

			return
				new ForceVector(vec, left.Unit);
		}

		/// <returns>
		///     A vector with components multiplied by a value
		/// </returns>
		public static ForceVector operator *(double multiplier, ForceVector vector) => new(vector.Values.Select(v => v * multiplier), vector.Unit);

		/// <inheritdoc cref="Matrix{T}.op_Multiply(Matrix{T}, Vector{T})"/>
		public static ForceVector operator *(Matrix<double> left, ForceVector right) => new(left * (Vector<double>) right, right.Unit);
		
		/// <inheritdoc cref="op_Multiply(double, ForceVector) " />
		public static ForceVector operator *(ForceVector vector, double multiplier) => multiplier * vector;

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static ForceVector operator -(ForceVector vector) => new(vector.Values.Select(v => -v), vector.Unit);


		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static ForceVector operator /(ForceVector vector, double divisor) => new(vector.Values.Select(v => v / divisor), vector.Unit);

		#endregion

	}
}