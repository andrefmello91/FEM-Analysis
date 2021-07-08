using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
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

		/// <inheritdoc cref="ICloneable.Clone" />
		public new ForceVector Clone() => (ForceVector) base.Clone();

		#endregion

		#region Operators

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Addition" />
		public static ForceVector operator +(ForceVector left, ForceVector right) => (ForceVector) ((ComponentVector<Force, ForceUnit>) left + right);

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Subtraction" />
		public static ForceVector operator -(ForceVector left, ForceVector right) => (ForceVector) ((ComponentVector<Force, ForceUnit>) left - right);

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Multiply(double,ComponentVector{TQuantity,TUnit}) " />
		public static ForceVector operator *(double multiplier, ForceVector vector) => (ForceVector) (multiplier * (ComponentVector<Force, ForceUnit>) vector);

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Multiply(double,ComponentVector{TQuantity,TUnit}) " />
		public static ForceVector operator *(ForceVector vector, double multiplier) => multiplier * vector;

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static ForceVector operator -(ForceVector vector) => (ForceVector) (-(ComponentVector<Force, ForceUnit>) vector);


		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static ForceVector operator /(ForceVector vector, double divisor) => (ForceVector) ((ComponentVector<Force, ForceUnit>) vector / divisor);

		#endregion

	}
}