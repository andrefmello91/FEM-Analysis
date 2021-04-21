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
		/// <inheritdoc/>
		public LinearAnalysis(FEMInput<IFiniteElement> femInput)
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

			// Assemble and simplify global stiffness and force vector
			UpdateStiffness();

			// Solve
			DisplacementVector = CalculateDisplacements(GlobalStiffness!, ForceVector)!;

			// Set displacements to grips
			FemInput.Grips.SetDisplacements(DisplacementVector);
			
			// Update element displacements
			FemInput.Elements.UpdateDisplacements();
			
			// Calculate element forces
			FemInput.Elements.CalculateForces();

			// Set Reactions
			FemInput.Grips.SetReactions(GetReactions());
		}
        #endregion

	}
}