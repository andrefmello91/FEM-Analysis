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
	public class ForceVector : ComponentVector<Force, ForceUnit>
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
			f.ConstraintIndex = femInput.ConstraintIndex;

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
			fi.ConstraintIndex = femInput.ConstraintIndex;

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
		public static ForceVector Zero(int size) => new(new double[size]);

#if NET5_0
		/// <inheritdoc />
		public override ForceVector Convert(ForceUnit unit) => new (Values.GetQuantities<Force, ForceUnit>(Unit).GetValues(unit), unit)
		{
			ConstraintIndex = ConstraintIndex
		};
		
		/// <inheritdoc cref="ICloneable.Clone" />
		public override ForceVector Clone() => new (Values, Unit)
		{
			ConstraintIndex = ConstraintIndex
		};

#else

		/// <inheritdoc />
		public override ComponentVector<Force, ForceUnit> Convert(ForceUnit unit) => new ForceVector(Values.GetQuantities<Force, ForceUnit>(Unit).GetValues(unit), unit)
		{
			ConstraintIndex = ConstraintIndex
		};

		/// <inheritdoc cref="ICloneable.Clone" />
		public override ComponentVector<Force, ForceUnit> Clone() => new ForceVector(Values, Unit)
		{
			ConstraintIndex = ConstraintIndex
		};

#endif

		/// <inheritdoc />
		public override Vector<double> Simplified() => Simplified(Tolerance);

		#endregion

		#region Operators

		/// <returns>
		///     A new vector with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ForceVector operator +(ForceVector left, ForceVector right) =>
			new(left.ToVector(left.Unit) + right.ToVector(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};

		/// <returns>
		///     A new vector with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentException">If left and right don't have the same dimensions.</exception>
		public static ForceVector operator -(ForceVector left, ForceVector right) =>
			new(left.ToVector(left.Unit) - right.ToVector(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};

		/// <returns>
		///     A vector with components multiplied by a value
		/// </returns>
		public static ForceVector operator *(double multiplier, ForceVector vector) =>
			new(vector.Values.Select(v => v * multiplier), vector.Unit)
			{
				ConstraintIndex = vector.ConstraintIndex
			};

		/// <inheritdoc cref="Matrix{T}.op_Multiply(Matrix{T}, Vector{T})"/>
		public static ForceVector operator *(Matrix<double> left, ForceVector right) =>
			new (left * right.ToVector(right.Unit), right.Unit)
			{
				ConstraintIndex = right.ConstraintIndex
			};
		
		/// <inheritdoc cref="op_Multiply(double, ForceVector) " />
		public static ForceVector operator *(ForceVector vector, double multiplier) => multiplier * vector;

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static ForceVector operator -(ForceVector vector) => new(vector.Select(v => -v.Value), vector.Unit)
		{
			ConstraintIndex = vector.ConstraintIndex
		};


		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static ForceVector operator /(ForceVector vector, double divisor) => new(vector.Select(v => v.Value / divisor), vector.Unit)
		{
			ConstraintIndex = vector.ConstraintIndex
		};

		#endregion

	}
}