using static andrefmello91.FEMAnalysis.StiffnessMatrix;

#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Linear analysis class.
	/// </summary>
	public class LinearAnalysis : Analysis<IFiniteElement>
	{

		#region Constructors

		/// <summary>
		///     Linear analysis constructor.
		/// </summary>
		/// <inheritdoc />
		public LinearAnalysis(IFEMInput<IFiniteElement> femInput)
			: base(femInput)
		{
		}

		#endregion

		#region Methods

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis{TFiniteElement}.ForceVector" /> (default: 1).</param>
		public void Execute(double loadFactor = 1)
		{
			// Set force vector
			ForceVector = FemInput.ForceVector * loadFactor;

			// Assemble stiffness
			UpdateStiffness();

			// Simplify global stiffness and force vector
			var stiffness = SimplifiedStiffness(GlobalStiffness!, FemInput.ConstraintIndex);
			var forces    = SimplifiedForces(ForceVector, FemInput.ConstraintIndex);

			// Solve
			DisplacementVector = CalculateDisplacements(stiffness, forces)!;

			// Set displacements to grips
			FemInput.Grips.SetDisplacements(DisplacementVector);

			// Update element displacements
			FemInput.UpdateDisplacements();

			// Calculate element forces
			FemInput.CalculateForces();

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}

		#endregion

	}
}