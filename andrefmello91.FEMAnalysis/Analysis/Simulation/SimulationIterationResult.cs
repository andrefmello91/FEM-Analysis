using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis.Simulation
{
	/// <summary>
	///		Iteration result class for simulations.
	/// </summary>
	internal class SimulationIterationResult : IterationResult, ICloneable<SimulationIterationResult>
	{
		/// <summary>
		///     The load factor increment of this iteration.
		/// </summary>
		public double LoadFactorIncrement { get; set; }

		/// <summary>
		///		The Arc-Length calculated for this iteration.
		/// </summary>
		public double ArcLength { get; set; }

		/// <summary>
		///     The displacement increment vector of this iteration.
		/// </summary>
		/// <returns>
		///		<see cref="IncrementFromResidual"/> + <see cref="LoadFactorIncrement"/> * <see cref="IncrementFromExternal"/>
		/// </returns>
		public new Vector<double> DisplacementIncrement => IncrementFromResidual + LoadFactorIncrement * IncrementFromExternal;

		/// <summary>
		///     The displacement increment vector from residual forces of this iteration.
		/// </summary>
		public Vector<double> IncrementFromResidual { get; set; }
		
		/// <summary>
		///     The displacement increment vector from external forces of this iteration.
		/// </summary>
		public Vector<double> IncrementFromExternal { get; set; }

		/// <summary>
		///		The parameter to set the sign of the load increments.
		/// </summary>
		public double StiffnessParameter { get; set; }
		
		/// <inheritdoc />
		internal SimulationIterationResult(int numberOfDoFs) : base(numberOfDoFs)
		{
		}

		/// <inheritdoc />
		internal SimulationIterationResult(Vector<double> displacements, Vector<double> residualForces, Matrix<double> stiffness) : base(displacements, residualForces, stiffness)
		{
		}

		/// <inheritdoc />
		public new SimulationIterationResult Clone() => new(Displacements.Clone(), ResidualForces.Clone(), Stiffness.Clone())
		{
			Number              = Number,
			LoadFactorIncrement = LoadFactorIncrement,
			ArcLength           = ArcLength,
			StiffnessParameter  = StiffnessParameter
		};
	}
}