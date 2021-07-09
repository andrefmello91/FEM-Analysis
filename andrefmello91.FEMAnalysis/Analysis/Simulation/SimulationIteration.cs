using andrefmello91.Extensions;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Iteration result class for simulations.
	/// </summary>
	public class SimulationIteration : Iteration, IIteration, ICloneable<SimulationIteration>
	{

		#region Properties

		/// <summary>
		///     The displacement increment vector from external forces of this iteration.
		/// </summary>
		public DisplacementVector IncrementFromExternal { get; private set; }

		/// <summary>
		///     The displacement increment vector from residual forces of this iteration.
		/// </summary>
		public DisplacementVector IncrementFromResidual { get; private set; }

		/// <summary>
		///     The load factor increment of this iteration.
		/// </summary>
		public double LoadFactorIncrement { get; set; }

		#region Interface Implementations

		/// <summary>
		///     The displacement increment vector of this iteration.
		/// </summary>
		/// <returns>
		///     <see cref="IncrementFromResidual" /> + <see cref="LoadFactorIncrement" /> * <see cref="IncrementFromExternal" />
		/// </returns>
		public override DisplacementVector DisplacementIncrement => IncrementFromResidual + LoadFactorIncrement * IncrementFromExternal;

		#endregion

		#endregion

		#region Constructors

		/// <inheritdoc />
		internal SimulationIteration(int numberOfDoFs)
			: base(numberOfDoFs)
		{
		}

		/// <inheritdoc />
		internal SimulationIteration(DisplacementVector displacements, ForceVector residualForces, StiffnessMatrix stiffness)
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
		/// <param name="updateDisplacements">Update displacement vector?</param>
		public void IncrementDisplacements(DisplacementVector? incrementFromResidual, DisplacementVector? incrementFromExternal, bool updateDisplacements = false)
		{
			if (incrementFromResidual is not null)
				IncrementFromResidual = incrementFromResidual;

			if (incrementFromExternal is not null)
				IncrementFromExternal = incrementFromExternal;

			if (updateDisplacements)
				UpdateDisplacements();
		}

		/// <summary>
		///     Add the displacement increment to displacement vector.
		/// </summary>
		public void UpdateDisplacements() => Displacements += DisplacementIncrement;

		#region Interface Implementations

		// /// <summary>
		// ///     Calculate the convergence of this iteration.
		// /// </summary>
		// /// <param name="initialIncrement">The displacement increment of the first iteration.</param>
		// public void CalculateConvergence(DisplacementVector initialIncrement) =>
		// 	DisplacementConvergence = NonlinearAnalysis.CalculateConvergence(DisplacementIncrement, initialIncrement);
		//
		// /// <inheritdoc/>
		// public override bool CheckConvergence(AnalysisParameters parameters) => base.CheckConvergence(parameters);

		/// <inheritdoc />
		IIteration ICloneable<IIteration>.Clone() => Clone();

		/// <inheritdoc />
		public new SimulationIteration Clone() => new((DisplacementVector) Displacements.Clone(), (ForceVector) ResidualForces.Clone(), (StiffnessMatrix) Stiffness.Clone())
		{
			Number                = Number,
			LoadFactorIncrement   = LoadFactorIncrement,
			IncrementFromResidual = IncrementFromResidual,
			IncrementFromExternal = IncrementFromExternal
		};

		#endregion

		#endregion

	}
}