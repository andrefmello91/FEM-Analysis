namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Secant simulation class.
	/// </summary>
	public class SecantSimulation : SecantAnalysis
	{

		/// <summary>
		///		Secant simulation constructor.
		/// </summary>
		/// <inheritdoc />
		public SecantSimulation(FEMInput femInput, double tolerance = 1E-06, int maxIterations = 10000, int minIterations = 2)
			: base(femInput, 50, tolerance, maxIterations, minIterations)
		{
		}

		/// <inheritdoc />
		protected override void StepAnalysis()
		{
			// Initiate first load step
			LoadStep = 1;

			while (true)
			{
				// Get the force vector
				CurrentForces = LoadFactor * ForceVector;

				// Iterate
				Iterate();

				// Verify if convergence was not reached
				if (Stop)
					return;

				// Set load step results
				SaveLoadStepResults();

				// Increment load step
				LoadStep++;
			}
		}
	}
}