using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Iteration result class for simulations.
	/// </summary>
	public class SimulationIteration : Iteration, IIteration, ICloneable<SimulationIteration>
	{

		#region Properties

		/// <summary>
		///     The Arc-Length calculated for this iteration.
		/// </summary>
		public double ArcLength { get; set; }

		/// <summary>
		///     The displacement increment vector of this iteration.
		/// </summary>
		/// <returns>
		///     <see cref="IncrementFromResidual" /> + <see cref="LoadFactorIncrement" /> * <see cref="IncrementFromExternal" />
		/// </returns>
		public override Vector<double> DisplacementIncrement => IncrementFromResidual + LoadFactorIncrement * IncrementFromExternal;

		/// <summary>
		///     The displacement increment vector from external forces of this iteration.
		/// </summary>
		public Vector<double> IncrementFromExternal { get; private set; }

		/// <summary>
		///     The displacement increment vector from residual forces of this iteration.
		/// </summary>
		public Vector<double> IncrementFromResidual { get; private set; }

		/// <summary>
		///     The load factor increment of this iteration.
		/// </summary>
		public double LoadFactorIncrement { get; set; }

		/// <summary>
		///     The parameter to set the sign of the load increments.
		/// </summary>
		public double StiffnessParameter { get; set; }

		#endregion

		#region Constructors

		/// <inheritdoc />
		internal SimulationIteration(int numberOfDoFs)
			: base(numberOfDoFs)
		{
		}

		/// <inheritdoc />
		internal SimulationIteration(Vector<double> displacements, Vector<double> residualForces, Matrix<double> stiffness)
			: base(displacements, residualForces, stiffness)
		{
		}

		#endregion

		#region Methods

		/// <summary>
		///     Increment displacements of this iteration.
		/// </summary>
		/// <param name="incrementFromResidual">The displacement increment vector from residual forces of this iteration.</param>
		/// <param name="incrementFromExternal">The displacement increment vector from external forces of this iteration.</param>
		public void IncrementDisplacements(Vector<double>? incrementFromResidual, Vector<double>? incrementFromExternal)
		{
			if (incrementFromResidual is not null)
				IncrementFromResidual = incrementFromResidual;

			if (incrementFromExternal is not null)
				IncrementFromExternal = incrementFromExternal;

			Displacements += DisplacementIncrement;
		}

		#region Interface Implementations

		/// <inheritdoc />
		public new SimulationIteration Clone() => new(Displacements.Clone(), ResidualForces.Clone(), Stiffness.Clone())
		{
			Number              = Number,
			LoadFactorIncrement = LoadFactorIncrement,
			ArcLength           = ArcLength,
			StiffnessParameter  = StiffnessParameter
		};
		
		/// <summary>
		///		Calculate the current stiffness parameter for defining the sign of the load factor increment.
		/// </summary>
		public void CalculateStiffnessParameter(LoadStep step)
		{
			if (Number <= 1)
			{
				StiffnessParameter = 1;
				return;
			}

			var inc = DisplacementIncrement;

			var k = (step.Forces.ToRowMatrix() * inc)[0] / (inc.ToRowMatrix() * inc)[0];

			StiffnessParameter = k / ((SimulationIteration) step.FirstIteration).StiffnessParameter;
		}

		/// <inheritdoc />
		IIteration ICloneable<IIteration>.Clone() => Clone();

		#endregion

		#endregion

	}
}