namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     The analysis types.
	/// </summary>
	public enum AnalysisType
	{
		/// <summary>
		///     Linear-elastic analysis.
		/// </summary>
		Linear,

		/// <summary>
		///     Nonlinear analysis.
		/// </summary>
		Nonlinear
	}

	/// <summary>
	///     Nonlinear solution procedures.
	/// </summary>
	public enum NonLinearSolver
	{
		/// <summary>
		///     Newton-Raphson nonlinear solver.
		/// </summary>
		NewtonRaphson,

		/// <summary>
		///     Modified Newton-Raphson nonlinear solver.
		/// </summary>
		ModifiedNewtonRaphson,

		/// <summary>
		///     Secant Method nonlinear solver.
		/// </summary>
		Secant
	}
}