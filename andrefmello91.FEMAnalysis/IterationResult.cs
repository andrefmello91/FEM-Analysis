using System.Collections.Generic;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Class for nonlinear iteration results.
	/// </summary>
	public class IterationResult : ICloneable<IterationResult>
	{

		#region Properties

		/// <summary>
		///     The convergence of this iteration.
		/// </summary>
		public double ForceConvergence { get; private set; }

		/// <summary>
		///     The displacement convergence of this iteration.
		/// </summary>
		public double DisplacementConvergence { get; private set; }

		/// <summary>
		///     The displacement vector of this iteration.
		/// </summary>
		public Vector<double> Displacements { get; set; }
		
		/// <summary>
		///     The displacement increment vector of this iteration.
		/// </summary>
		public Vector<double> DisplacementIncrement { get; set; }

		/// <summary>
		///     The number of this iteration.
		/// </summary>
		public int Number { get; set; }

		/// <summary>
		///     The internal force vector of this iteration.
		/// </summary>
		public Vector<double> InternalForces { get; set; }
		
		/// <summary>
		///     The residual force vector of this iteration.
		/// </summary>
		public Vector<double> ResidualForces { get; private set; }

		/// <summary>
		///     The stiffness matrix of this iteration.
		/// </summary>
		public Matrix<double> Stiffness { get; set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Create an iteration object with elements composed by zero..
		/// </summary>
		/// <param name="numberOfDoFs">The number of degrees of freedom.</param>
		public IterationResult(int numberOfDoFs)
			: this(Vector<double>.Build.Dense(numberOfDoFs), Vector<double>.Build.Dense(numberOfDoFs), Matrix<double>.Build.Dense(numberOfDoFs, numberOfDoFs))
		{
		}

		/// <summary>
		///     Create an iteration object.
		/// </summary>
		/// <param name="displacements">The displacement vector of this iteration.</param>
		/// <param name="residualForces">The residual force vector of this iteration.</param>
		/// <param name="stiffness">The stiffness matrix of this iteration.</param>
		public IterationResult(Vector<double> displacements, Vector<double> residualForces, Matrix<double> stiffness)
		{
			Displacements  = displacements;
			ResidualForces = residualForces;
			Stiffness      = stiffness;
			InternalForces = Vector<double>.Build.Dense(displacements.Count);
		}

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the convergence of this iteration.
		/// </summary>
		/// <param name="appliedForces">The applied forces of the current step.</param>
		public void CalculateForceConvergence(IEnumerable<double> appliedForces) =>
			ForceConvergence = NonlinearAnalysis.CalculateConvergence(ResidualForces, appliedForces);

		/// <summary>
		///     Calculate the displacement convergence of this iteration.
		/// </summary>
		/// <param name="initialIncrement">The displacement increment of the first iteration.</param>
		public void CalculateDisplacementConvergence(IEnumerable<double> initialIncrement) =>
			DisplacementConvergence = NonlinearAnalysis.CalculateConvergence(DisplacementIncrement, initialIncrement);

		/// <summary>
		///		Update forces in this iteration.
		/// </summary>
		/// <param name="appliedForces">The vector of applied forces of the current step.</param>
		/// <param name="internalForces">The vector of internal forces.</param>
		public void UpdateForces(Vector<double> appliedForces, Vector<double> internalForces)
		{
			InternalForces = internalForces;
			ResidualForces = internalForces - appliedForces;
		}
		
		#region Interface Implementations

		/// <inheritdoc />
		public IterationResult Clone() => new(Displacements.Clone(), ResidualForces.Clone(), Stiffness.Clone()) { Number = Number };

		#endregion

		#endregion

		#region Operators

		/// <summary>
		///     Get the number of a iteration.
		/// </summary>
		/// <returns>
		///     <see cref="IterationResult.Number" />
		/// </returns>
		public static explicit operator int(IterationResult iterationResult) => iterationResult.Number;

		#endregion

	}
}