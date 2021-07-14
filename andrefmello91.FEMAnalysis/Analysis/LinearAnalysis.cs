using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Linear analysis class.
	/// </summary>
	public class LinearAnalysis : Analysis
	{

		#region Constructors

		/// <summary>
		///     Linear analysis constructor.
		/// </summary>
		/// <inheritdoc />
		public LinearAnalysis(IFEMInput femInput)
			: base(femInput)
		{
		}

		#endregion

		#region Methods

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.Forces" /> (default: 1).</param>
		public void Execute(double loadFactor = 1)
		{
			// Set force vector
			if (!loadFactor.Approx(1))
				Forces = (ForceVector) (Forces * loadFactor);

			// Solve
			Displacements = GlobalStiffness.Solve(Forces);

			// Set displacements to grips
			FemInput.Grips.SetDisplacements(Displacements);

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