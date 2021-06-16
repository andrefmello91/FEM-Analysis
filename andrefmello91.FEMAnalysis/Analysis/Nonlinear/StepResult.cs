using System.Collections.Generic;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Class for step results.
	/// </summary>
	public class StepResult : List<IterationResult>, ICloneable<StepResult>
	{

		#region Properties

		/// <summary>
		///     The convergence of this step.
		/// </summary>
		public double Convergence { get; set; }
		
		/// <summary>
		///     The load factor of this step.
		/// </summary>
		public double LoadFactor { get; set; }

		/// <summary>
		///     The displacement vector of this step.
		/// </summary>
		public Vector<double> Displacements { get; set; }

		/// <summary>
		///     The force vector of this step.
		/// </summary>
		public Vector<double> Forces { get; set; }

		/// <summary>
		///     The status of this step. True if it was calculated.
		/// </summary>
		public bool IsCalculated { get; set; }

		/// <summary>
		///     The monitored displacement of this step.
		/// </summary>
		public MonitoredDisplacement? MonitoredDisplacement { get; set; }

		/// <summary>
		///     The number of this step.
		/// </summary>
		public int Number { get; set; }

		/// <summary>
		///     The stiffness matrix of this step.
		/// </summary>
		public Matrix<double> Stiffness { get; set; }
		
		#endregion

		#region Constructors

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="number">The number of this step.</param>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		public StepResult(int numberOfDoFs, int number = 0)
			: this(number, Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs))
		{
		}

		/// <summary>
		///     Create a step object.
		/// </summary>
		/// <param name="number">The number of this step.</param>
		/// <param name="forces">The force vector of this step.</param>
		public StepResult(Vector<double> forces, int number = 0)
			: this(number, forces, Vector<double>.Build.Dense(forces.Count), Matrix<double>.Build.Dense(forces.Count, forces.Count))
		{
		}

		/// <inheritdoc cref="StepResult" />
		/// <param name="displacements">The displacement vector of this step.</param>
		/// <param name="stiffness">The stiffness matrix of this step.</param>
		public StepResult(int number, Vector<double> forces, Vector<double> displacements, Matrix<double> stiffness)
		{
			Number        = number;
			Forces        = forces;
			Displacements = displacements;
			Stiffness     = stiffness;
		}

		#endregion

		#region Methods

		#region Interface Implementations

		/// <inheritdoc />
		public StepResult Clone() => new(Number, Forces.Clone(), Displacements.Clone(), Stiffness.Clone())
		{
			LoadFactor = LoadFactor
		};

		/// <inheritdoc />
		public override string ToString() => $"Step {Number}";

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the number of a step.
		/// </summary>
		/// <returns>
		///     <see cref="StepResult.Number" />
		/// </returns>
		public static explicit operator int(StepResult stepResult) => stepResult.Number;

		/// <summary>
		///		Check the step number.
		/// </summary>
		/// <returns>
		///		True if the step number is equal to the right number.
		/// </returns>
		public static bool operator ==(StepResult left, int right) => left.Number == right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is not equal to the right number.
		/// </returns>
		public static bool operator !=(StepResult left, int right) => left.Number != right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is bigger than the right number.
		/// </returns>
		public static bool operator >(StepResult left, int right) => left.Number > right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is smaller than the right number.
		/// </returns>
		public static bool operator <(StepResult left, int right) => left.Number < right;
		
		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is bigger or equal to the right number.
		/// </returns>
		public static bool operator >=(StepResult left, int right) => left.Number >= right;

		/// <inheritdoc cref="op_Equality"/>
		/// <returns>
		///		True if the step number is smaller or equal to the right number.
		/// </returns>
		public static bool operator <=(StepResult left, int right) => left.Number <= right;

		#endregion

	}
}