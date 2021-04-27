using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Class for load step results.
	/// </summary>
	public class LoadStepResult : ICloneable<LoadStepResult>
	{

		#region Properties

		/// <summary>
		///     The convergence of this load step.
		/// </summary>
		public double Convergence { get; set; }

		/// <summary>
		///     The displacement vector of this load step.
		/// </summary>
		public Vector<double> Displacements { get; set; }

		/// <summary>
		///     The force vector of this load step.
		/// </summary>
		public Vector<double> Forces { get; set; }

		/// <summary>
		///     The status of this load step. True if it was calculated.
		/// </summary>
		public bool IsCalculated { get; set; }

		/// <summary>
		///     The monitored displacement of this load step.
		/// </summary>
		public MonitoredDisplacement? MonitoredDisplacement { get; set; }

		/// <summary>
		///     The number of this load step.
		/// </summary>
		public int Number { get; set; }

		/// <summary>
		///     The stiffness matrix of this load step.
		/// </summary>
		public Matrix<double> Stiffness { get; set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Create a load step object.
		/// </summary>
		/// <param name="number">The number of this load step.</param>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		public LoadStepResult(int numberOfDoFs, int number = 0)
			: this(number, Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs))
		{
		}

		/// <summary>
		///     Create a load step object.
		/// </summary>
		/// <param name="number">The number of this load step.</param>
		/// <param name="forces">The force vector of this load step.</param>
		public LoadStepResult(Vector<double> forces, int number = 0)
			: this(number, forces, Vector<double>.Build.Dense(forces.Count), Matrix<double>.Build.Dense(forces.Count, forces.Count))
		{
		}

		/// <inheritdoc cref="LoadStepResult(int, Vector{double})" />
		/// <param name="displacements">The displacement vector of this load step.</param>
		/// <param name="stiffness">The stiffness matrix of this load step.</param>
		public LoadStepResult(int number, Vector<double> forces, Vector<double> displacements, Matrix<double> stiffness)
		{
			Number        = number;
			Forces        = forces;
			Displacements = displacements;
			Stiffness     = stiffness;
		}

		#endregion

		#region Methods

		/// <inheritdoc />
		public LoadStepResult Clone() => new(Number, Forces.Clone(), Displacements.Clone(), Stiffness.Clone());

		#endregion

		#region Operators

		/// <summary>
		///     Get the number of a load step.
		/// </summary>
		/// <returns>
		///     <see cref="LoadStepResult.Number" />
		/// </returns>
		public static explicit operator int(LoadStepResult loadStepResult) => loadStepResult.Number;

		#endregion

	}
}